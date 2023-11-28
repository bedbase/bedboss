from typing import Union
import json
import yaml
import os
import requests
import pypiper
import bbconf
import logging
from geniml.io import RegionSet


_LOGGER = logging.getLogger("bedboss")

SCHEMA_PATH_BEDSTAT = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "pep_schema.yaml"
)


def convert_unit(size_in_bytes: int) -> str:
    """
    Convert the size from bytes to other units like KB, MB or GB
    :param int size_in_bytes: size in bytes
    :return str: File size as string in different units
    """
    if size_in_bytes < 1024:
        return str(size_in_bytes) + "bytes"
    elif size_in_bytes in range(1024, 1024 * 1024):
        return str(round(size_in_bytes / 1024, 2)) + "KB"
    elif size_in_bytes in range(1024 * 1024, 1024 * 1024 * 1024):
        return str(round(size_in_bytes / (1024 * 1024))) + "MB"
    elif size_in_bytes >= 1024 * 1024 * 1024:
        return str(round(size_in_bytes / (1024 * 1024 * 1024))) + "GB"


def bedstat(
    bedfile: str,
    bedbase_config: Union[str, bbconf.BedBaseConf],
    genome: str,
    outfolder: str,
    ensdb: str = None,
    open_signal_matrix: str = None,
    bigbed: str = None,
    treatment: str = None,
    description: str = None,
    cell_type: str = None,
    other_metadata: dict = None,
    just_db_commit: bool = False,
    no_db_commit: bool = False,
    force_overwrite: bool = False,
    skip_qdrant: bool = True,
    pm: pypiper.PipelineManager = None,
    **kwargs,
) -> str:
    """
    Run bedstat pipeline - pipeline for obtaining statistics about bed files
        and inserting them into the database

    :param str bedfile: the full path to the bed file to process
    :param str bigbed: the full path to the bigbed file. Defaults to None.
        (bigbed won't be created and some producing of some statistics will
        be skipped.)
    :param str bedbase_config: The path to the bedbase configuration file, or bbconf object
    :param str open_signal_matrix: a full path to the openSignalMatrix
        required for the tissue specificity plots
    :param str outfolder: The folder for storing the pipeline results.
    :param str genome: genome assembly of the sample
    :param str ensdb: a full path to the ensdb gtf file required for genomes
        not in GDdata
    :param str description: a description of the bed file
    :param str treatment: a treatment of the bed file
    :param str cell_type: a cell type of the bed file
    :param dict other_metadata: a dictionary of other metadata to pass
    :param bool just_db_commit: whether just to commit the JSON to the database
    :param bool no_db_commit: whether the JSON commit to the database should be
        skipped
    :param skip_qdrant: whether to skip qdrant indexing [Default: True]
    :param bool force_overwrite: whether to overwrite the existing record
    :param pm: pypiper object

    :return: bed_digest: the digest of the bed file
    """
    # TODO why are we no longer using bbconf to get the output path?
    # outfolder_stats = bbc.get_bedstat_output_path()
    outfolder_stats = os.path.join(outfolder, "output", "bedstat_output")
    try:
        os.makedirs(outfolder_stats)
    except FileExistsError:
        pass

    # if bbconf is a string, create a bbconf object
    if isinstance(bedbase_config, str):
        bbc = bbconf.BedBaseConf(config_path=bedbase_config, database_only=True)
    else:
        bbc = bedbase_config

    bed_digest = RegionSet(bedfile).identifier
    bedfile_name = os.path.split(bedfile)[1]

    fileid = os.path.splitext(os.path.splitext(bedfile_name)[0])[0]
    outfolder = os.path.abspath(os.path.join(outfolder_stats, bed_digest))
    try:
        os.makedirs(outfolder)
    except FileExistsError:
        pass
    json_file_path = os.path.abspath(os.path.join(outfolder, fileid + ".json"))
    json_plots_file_path = os.path.abspath(
        os.path.join(outfolder, fileid + "_plots.json")
    )
    bed_relpath = os.path.relpath(
        bedfile,
        os.path.abspath(os.path.join(outfolder_stats, os.pardir, os.pardir)),
    )
    bigbed_relpath = os.path.relpath(
        os.path.join(bigbed, fileid + ".bigBed"),
        os.path.abspath(os.path.join(outfolder_stats, os.pardir, os.pardir)),
    )
    if not just_db_commit:
        if not pm:
            pm_out_path = os.path.abspath(
                os.path.join(outfolder_stats, "pypiper", bed_digest)
            )
            try:
                os.makedirs(pm_out_path)
            except FileExistsError:
                pass
            pm = pypiper.PipelineManager(
                name="bedstat-pipeline",
                outfolder=pm_out_path,
                pipestat_sample_name=bed_digest,
            )
            stop_pipeline = True
        else:
            stop_pipeline = False

        rscript_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "bedstat",
            "tools",
            "regionstat.R",
        )
        assert os.path.exists(rscript_path), FileNotFoundError(
            f"'{rscript_path}' script not found"
        )
        command = (
            f"Rscript {rscript_path} --bedfilePath={bedfile} "
            f"--fileId={fileid} --openSignalMatrix={open_signal_matrix} "
            f"--outputFolder={outfolder} --genome={genome} "
            f"--ensdb={ensdb} --digest={bed_digest}"
        )

        pm.run(cmd=command, target=json_file_path)

    # commit to the database if no_db_commit is not set
    if not no_db_commit:
        data = {}
        if os.path.exists(json_file_path):
            with open(json_file_path, "r", encoding="utf-8") as f:
                data = json.loads(f.read())
        if os.path.exists(json_plots_file_path):
            with open(json_plots_file_path, "r", encoding="utf-8") as f_plots:
                plots = json.loads(f_plots.read())
        else:
            plots = []

        if not other_metadata:
            other_metadata = {}
        other_metadata.update(
            {
                "description": description,
                "treatment": treatment,
                "cell_type": cell_type,
            }
        )

        # unlist the data, since the output of regionstat.R is a dict of lists of
        # length 1 and force keys to lower to correspond with the
        # postgres column identifiers
        data = {k.lower(): v[0] if isinstance(v, list) else v for k, v in data.items()}
        data.update(
            {
                "bedfile": {
                    "path": bed_relpath,
                    "size": convert_unit(os.path.getsize(bedfile)),
                    "title": "Path to the BED file",
                }
            }
        )

        if os.path.exists(os.path.join(bigbed, fileid + ".bigBed")):
            data.update(
                {
                    "bigbedfile": {
                        "path": bigbed_relpath,
                        "size": convert_unit(
                            os.path.getsize(os.path.join(bigbed, fileid + ".bigBed"))
                        ),
                        "title": "Path to the big BED file",
                    }
                }
            )

            if not os.path.islink(os.path.join(bigbed, fileid + ".bigBed")):
                digest = requests.get(
                    f"http://refgenomes.databio.org/genomes/genome_digest/{genome}"
                ).text.strip('""')

                data.update(
                    {
                        "genome": {
                            "alias": genome,
                            "digest": digest,
                        }
                    }
                )
        else:
            data.update(
                {
                    "genome": {
                        "alias": genome,
                        "digest": "",
                    }
                }
            )

        for plot in plots:
            plot_id = plot["name"]
            del plot["name"]
            data.update({plot_id: plot})

        # deleting md5sum, because it is record_identifier
        del data["md5sum"]

        # add added_to_qdrant to the data
        data["other"] = other_metadata
        data["added_to_qdrant"] = False

        bbc.bed.report(
            record_identifier=bed_digest,
            values=data,
            force_overwrite=force_overwrite,
        )

    if not skip_qdrant:
        bbc.add_bed_to_qdrant(
            bed_id=bed_digest,
            bed_file=bedfile,
            payload={"fileid": fileid},
        )
        bbc.bed.report(
            record_identifier=bed_digest,
            values={"added_to_qdrant": True},
            force_overwrite=True,
        )

    if stop_pipeline:
        pm.stop_pipeline()
    return bed_digest
