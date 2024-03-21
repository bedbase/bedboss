import logging
import os
from typing import NoReturn, Union

import pypiper
from argparse import Namespace
import logmuse
import peppy
from eido import validate_project
import bbconf
import subprocess

import pephubclient
from pephubclient.helpers import is_registry_path
from bbconf.bbagent import BedBaseAgent


from bedboss.bedstat.bedstat import bedstat
from bedboss.bedmaker.bedmaker import make_all
from bedboss.bedqc.bedqc import bedqc
from bedboss.bedbuncher import run_bedbuncher
from bedboss.qdrant_index import add_to_qdrant
from bedboss.cli import build_argparser
from bedboss.const import (
    BEDBOSS_PEP_SCHEMA_PATH,
)
from bedboss.models import (
    BedStatCLIModel,
    BedMakerCLIModel,
    BedQCCLIModel,
    StatsUpload,
    PlotsUpload,
    FilesUpload,
    BedClassificationUpload,
)
from bedboss.utils import (
    standardize_genome_name,
    get_genome_digest,
)
from bedboss.exceptions import BedBossException
from bedboss._version import __version__

_LOGGER = logging.getLogger("bedboss")


def requirements_check() -> None:
    """
    Check if all requirements are installed

    :return: None
    """
    _LOGGER.info("Checking requirements...")
    subprocess.run(
        ["bash", f"{os.path.dirname(os.path.abspath(__file__))}/requirements_test.sh"]
    )


def run_all(
    input_file: str,
    input_type: str,
    outfolder: str,
    genome: str,
    bedbase_config: str,
    rfg_config: str = None,
    narrowpeak: bool = False,
    check_qc: bool = True,
    chrom_sizes: str = None,
    open_signal_matrix: str = None,
    ensdb: str = None,
    other_metadata: dict = None,
    just_db_commit: bool = False,
    db_commit: bool = True,
    force_overwrite: bool = False,
    upload_qdrant: bool = False,
    upload_s3: bool = False,
    upload_pephub: bool = False,
    pm: pypiper.PipelineManager = None,
) -> str:
    """
    Run bedboss: bedmaker -> bedqc -> bedclassifier -> bedstat -> upload to s3, qdrant, pephub, and bedbase.

    :param str input_file: Input file [required]
    :param str input_type: Input type [required] options: (bigwig|bedgraph|bed|bigbed|wig)
    :param str outfolder: Folder, where output should be saved  [required]
    :param str genome: genome_assembly of the sample. [required] options: (hg19, hg38, mm10) # TODO: add more
    :param Union[str, bbconf.BedBaseConf] bedbase_config: The path to the bedbase configuration file, or bbconf object.
    :param str rfg_config: file path to the genome config file [optional]
    :param bool narrowpeak: whether the regions are narrow
        (transcription factor implies narrow, histone mark implies broad peaks) [optional]
    :param bool check_qc: set True to run quality control during badmaking [optional] (default: True)
    :param str chrom_sizes: a full path to the chrom.sizes required for the bedtobigbed conversion [optional]
    :param str open_signal_matrix: a full path to the openSignalMatrix required for the tissue [optional]
    :param dict other_metadata: a dict containing all attributes from the sample
    :param str ensdb: a full path to the ensdb gtf file required for genomes not in GDdata [optional]
        (basically genomes that's not in GDdata)
    :param bool just_db_commit: whether just to commit the JSON to the database (default: False)
    :param bool force_overwrite: force overwrite analysis

    :param bool db_commit: whether the JSON commit to the database should be skipped (default: False)
    :param bool upload_qdrant: whether to skip qdrant indexing
    :param bool upload_s3: whether to upload to s3
    :param bool upload_pephub: whether to push bedfiles and metadata to pephub (default: False)
    :param pypiper.PipelineManager pm: pypiper object
    :return str bed_digest: bed digest
    """
    bbagent = BedBaseAgent(bedbase_config)

    genome = standardize_genome_name(genome)

    _LOGGER.info(f"Input file = '{input_file}'")
    _LOGGER.info(f"Output folder = '{outfolder}'")

    if not pm:
        pm_out_folder = os.path.join(os.path.abspath(outfolder), "pipeline_manager")
        _LOGGER.info(f"Pipeline info folder = '{pm_out_folder}'")
        pm = pypiper.PipelineManager(
            name="bedboss-pipeline",
            outfolder=pm_out_folder,
            version=__version__,
            recover=True,
        )
        stop_pipeline = True
    else:
        stop_pipeline = False

    bed_metadata = make_all(
        input_file=input_file,
        input_type=input_type,
        output_path=outfolder,
        genome=genome,
        rfg_config=rfg_config,
        narrowpeak=narrowpeak,
        check_qc=check_qc,
        chrom_sizes=chrom_sizes,
        pm=pm,
    )
    if not other_metadata:
        other_metadata = {}

    statistics_dict = bedstat(
        bedfile=bed_metadata.bed_file,
        outfolder=outfolder,
        genome=genome,
        ensdb=ensdb,
        bed_digest=bed_metadata.bed_digest,
        open_signal_matrix=open_signal_matrix,
        just_db_commit=just_db_commit,
        pm=pm,
    )
    statistics_dict["bed_type"] = bed_metadata.bed_type
    statistics_dict["bed_format"] = bed_metadata.bed_format.value

    if bed_metadata.bigbed_file:
        genome_digest = get_genome_digest(genome)
    else:
        genome_digest = None

    stats = StatsUpload(**statistics_dict)
    plots = PlotsUpload(**statistics_dict)

    files = FilesUpload(
        bedfile={
            "name": "bedfile",
            "path": bed_metadata.bed_file,
            "description": "Path to the BED file",
        },
        bigbedfile={
            "name": "bigbedfile",
            "path": bed_metadata.bigbed_file,
            "description": "Path to the bigbed file",
        },
    )

    classification = BedClassificationUpload(
        name=bed_metadata.bed_digest,
        genome_digest=genome_digest,
        genome_alias=genome,
        bed_type=bed_metadata.bed_type,
        bed_format=bed_metadata.bed_format.value,
    )

    bbagent.bed.add(
        identifier=bed_metadata.bed_digest,
        stats=stats.model_dump(exclude_unset=True),
        metadata=other_metadata,
        plots=plots.model_dump(exclude_unset=True),
        files=files.model_dump(exclude_unset=True),
        classification=classification.model_dump(exclude_unset=True),
        add_to_qdrant=upload_qdrant,
        upload_pephub=upload_pephub,
        upload_s3=upload_s3,
        local_path=outfolder,
        overwrite=force_overwrite,
        nofail=True
    )

    if stop_pipeline:
        pm.stop_pipeline()

    return bed_metadata.bed_digest


def insert_pep(
    bedbase_config: str,
    output_folder: str,
    pep: Union[str, peppy.Project],
    rfg_config: str = None,
    create_bedset: bool = True,
    check_qc: bool = True,
    ensdb: str = None,
    db_commit: bool = True,
    just_db_commit: bool = False,
    force_overwrite: bool = False,
    upload_s3: bool = False,
    upload_pephub: bool = False,
    upload_qdrant: bool = False,
    pm: pypiper.PipelineManager = None,
    *args,
    **kwargs,
) -> None:
    """
    Run all bedboss pipelines for all samples in the pep file.
    bedmaker -> bedqc -> bedstat -> qdrant_indexing -> bedbuncher

    :param str bedbase_config: bedbase configuration file path
    :param str output_folder: output statistics folder
    :param Union[str, peppy.Project] pep: path to the pep file or pephub registry path
    :param str rfg_config: path to the genome config file (refgenie)
    :param bool create_bedset: whether to create bedset
    :param bool upload_qdrant: whether to upload bedfiles to qdrant
    :param bool check_qc: whether to run quality control during badmaking
    :param str ensdb: a full path to the ensdb gtf file required for genomes not in GDdata
    :param bool just_db_commit: whether save only to the database (Without saving locally )
    :param bool db_commit: whether to upload data to the database
    :param bool force_overwrite: whether to overwrite the existing record
    :param bool upload_s3: whether to upload to s3
    :param bool upload_pephub: whether to push bedfiles and metadata to pephub (default: False)
    :param bool upload_qdrant: whether to execute qdrant indexing
    :param pypiper.PipelineManager pm: pypiper object
    :return: None
    """

    _LOGGER.warning(f"!Unused arguments: {kwargs}")
    failed_samples = []
    pephub_registry_path = None
    if isinstance(pep, peppy.Project):
        pass
    elif isinstance(pep, str):
        if is_registry_path(pep):
            pephub_registry_path = pep
            pep = pephubclient.PEPHubClient().load_project(pep)
        else:
            pep = peppy.Project(pep)
    else:
        raise BedBossException("Incorrect pep type. Exiting...")

    bbc = bbconf.BedBaseConf(config_path=bedbase_config, database_only=True)

    validate_project(pep, BEDBOSS_PEP_SCHEMA_PATH)

    for i, pep_sample in enumerate(pep.samples):
        _LOGGER.info(f"Running bedboss pipeline for {pep_sample.sample_name}")
        if pep_sample.get("file_type"):
            if pep_sample.get("file_type").lower() == "narrowpeak":
                is_narrow_peak = True
            else:
                is_narrow_peak = False
        else:
            is_narrow_peak = False
        try:
            bed_id = run_all(
                input_file=pep_sample.input_file,
                input_type=pep_sample.input_type,
                genome=pep_sample.genome,
                narrowpeak=is_narrow_peak,
                chrom_sizes=pep_sample.get("chrom_sizes"),
                open_signal_matrix=pep_sample.get("open_signal_matrix"),
                other_metadata=pep_sample.to_dict(),
                outfolder=output_folder,
                bedbase_config=bbc,
                rfg_config=rfg_config,
                check_qc=check_qc,
                ensdb=ensdb,
                just_db_commit=just_db_commit,
                db_commit=db_commit,
                force_overwrite=force_overwrite,
                upload_qdrant=upload_qdrant,
                upload_s3=upload_s3,
                upload_pephub=upload_pephub,
                pm=pm,
            )
            pep.samples[i].record_identifier = bed_id
        except BedBossException as e:
            _LOGGER.error(f"Failed to process {pep_sample.sample_name}. See {e}")
            failed_samples.append(pep_sample.sample_name)

    else:
        _LOGGER.info("Skipping uploading to s3. Flag `upload_s3` is set to False")

    if create_bedset:
        _LOGGER.info(f"Creating bedset from {pep.name}")
        run_bedbuncher(
            bedbase_config=bbc,
            bedset_pep=pep,
            pephub_registry_path=pephub_registry_path,
            upload_pephub=upload_pephub,
        )
    else:
        _LOGGER.info(
            f"Skipping bedset creation. Create_bedset is set to {create_bedset}"
        )
    _LOGGER.info(f"Failed samples: {failed_samples}")


def main(test_args: dict = None) -> NoReturn:
    """
    Run pipeline that was specified in as positional argument.

    :param str test_args: one of the bedboss pipelines
    """
    parser = build_argparser()
    if test_args:
        args = Namespace(**test_args)
    else:
        args, _ = parser.parse_known_args()
        global _LOGGER
        _LOGGER = logmuse.logger_via_cli(args, make_root=True)

    args_dict = vars(args)

    pm_out_folder = (
        args_dict.get("outfolder")
        or args_dict.get("output_folder")
        or "test_outfolder",
    )
    pm_out_folder = os.path.join(os.path.abspath(pm_out_folder[0]), "pipeline_manager")
    pm = pypiper.PipelineManager(
        name="bedboss-pipeline",
        outfolder=pm_out_folder,
        version=__version__,
        # args=args,
        multi=args_dict.get("multy", False),
        recover=True,
    )
    if args_dict["command"] == "all":
        run_all(pm=pm, **args_dict)
    elif args_dict["command"] == "insert":
        insert_pep(pm=pm, **args_dict)
    elif args_dict["command"] == "make":
        make_all(**BedMakerCLIModel(pm=pm, **args_dict).model_dump())
    elif args_dict["command"] == "qc":
        bedqc(**BedQCCLIModel(pm=pm, **args_dict).model_dump())
    elif args_dict["command"] == "stat":
        bedstat(**BedStatCLIModel(pm=pm, **args_dict).model_dump())
    elif args_dict["command"] == "bunch":
        run_bedbuncher(pm=pm, **args_dict)
    elif args_dict["command"] == "index":
        add_to_qdrant(pm=pm, **args_dict)
    elif args_dict["command"] == "requirements-check":
        requirements_check()
    else:
        parser.print_help()
        # raise Exception("Incorrect pipeline name.")
    pm.stop_pipeline()
