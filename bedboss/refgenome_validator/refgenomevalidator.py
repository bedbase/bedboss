from typing import Optional, Union, List, Dict
import os
import subprocess

from bedboss.exceptions import ValidatorException
from bedboss.refgenome_validator.genome_model import GenomeModel
from bedboss.refgenome_validator.models import (
    ChromNameStats,
    ChromLengthStats,
    SequenceFitStats,
    CompatibilityStats,
    RatingModel,
    CompatibilityConcise,
)
from bedboss.refgenome_validator.utils import get_bed_chrom_info

try:
    IGD_LOCATION = os.environ["IGD_LOCATION"]
except:
    # Local installation of C version of IGD
    IGD_LOCATION = f"/home/drc/GITHUB/igd/IGD/bin/igd"


class ReferenceValidator:
    """
    This is primary class for creating a compatibility vector
    An object of this class is to be created once and then used for the entirety of a pipeline,e.g. Bedboss.
    """

    def __init__(
        self,
        genome_models: Optional[List[GenomeModel]] = None,
        igd_path: Optional[str] = None,
    ):
        """
        Initialization method

        :param genome_models: this is a list of GenomeModels that will be checked against a bed file
        :param igd_path: path to a local IGD file containing ALL excluded ranges intervals for IGD overlap assessment, if not provided these metrics are not computed.
        """

        if not genome_models:
            genome_models = self._build_default_models()
        elif isinstance(genome_models, str):
            genome_models = list(genome_models)
        elif not isinstance(genome_models, list):
            raise ValidatorException(
                reason="A list of GenomeModels must be provided to initialize the Validator class"
            )

        self.genome_models: List[GenomeModel] = genome_models
        self.igd_path = igd_path

    @staticmethod
    def calculate_chrom_stats(
        bed_chrom_sizes: dict, genome_chrom_sizes: dict
    ) -> CompatibilityStats:
        """
        Given two dicts of chroms (key) and their sizes (values)
        determine overlap and sequence fit

        Calculates Stats associated with comparison of chrom names, chrom lengths, and sequence fits

        :param dict bed_chrom_sizes: dict of a bedfile's chrom size
        :param dict genome_chrom_sizes: dict of a GenomeModel's chrom sizes

        return dict: returns a dictionary with information on Query vs Model, e.g. chrom names QueryvsModel
        """

        # Layer 1: Check names and Determine XS (Extra Sequences) via Calculation of Recall/Sensitivity
        # Q = Query, M = Model
        q_and_m = 0  # how many chrom names are in the query and the genome model?
        q_and_not_m = (
            0  # how many chrom names are in the query but not in the genome model?
        )
        not_q_and_m = (
            0  # how many chrom names are not in the query but are in the genome model?
        )
        passed_chrom_names = True  # does this bed file pass the first layer of testing?

        query_keys_present = []  # These keys are used for seq fit calculation
        for key in list(bed_chrom_sizes.keys()):
            if key not in genome_chrom_sizes:
                q_and_not_m += 1
            if key in genome_chrom_sizes:
                q_and_m += 1
                query_keys_present.append(key)
        for key in list(genome_chrom_sizes.keys()):
            if key not in bed_chrom_sizes:
                not_q_and_m += 1

        # Calculate the Jaccard Index for Chrom Names
        bed_chrom_set = set(list(bed_chrom_sizes.keys()))
        genome_chrom_set = set(list(genome_chrom_sizes.keys()))
        chrom_intersection = bed_chrom_set.intersection(genome_chrom_set)
        chrom_union = bed_chrom_set.union(chrom_intersection)
        chrom_jaccard_index = len(chrom_intersection) / len(chrom_union)

        # Alternative Method for Calculating Jaccard_index for binary classification
        # JI = TP/(TP+FP+FN)
        jaccard_binary = q_and_m / (q_and_m + not_q_and_m + q_and_not_m)

        # What is our threshold for passing layer 1?
        if q_and_not_m > 0:
            passed_chrom_names = False

        # Calculate sensitivity for chrom names
        # defined as XS -> Extra Sequences
        sensitivity = q_and_m / (q_and_m + q_and_not_m)

        # Assign Stats
        name_stats = ChromNameStats(
            xs=sensitivity,
            q_and_m=q_and_m,
            q_and_not_m=q_and_not_m,
            not_q_and_m=not_q_and_m,
            jaccard_index=chrom_jaccard_index,
            jaccard_index_binary=jaccard_binary,
            passed_chrom_names=passed_chrom_names,
        )

        # Layer 2:  Check Lengths, but only if layer 1 is passing
        if passed_chrom_names:
            chroms_beyond_range = False
            num_of_chrom_beyond = 0
            num_chrom_within_bounds = 0

            for key in list(bed_chrom_sizes.keys()):
                if key in genome_chrom_sizes:
                    if bed_chrom_sizes[key] > genome_chrom_sizes[key]:
                        num_of_chrom_beyond += 1
                        chroms_beyond_range = True
                    else:
                        num_chrom_within_bounds += 1

            # Calculate recall/sensitivity for chrom lengths
            # defined as OOBR -> Out of Bounds Range
            sensitivity = num_chrom_within_bounds / (
                num_chrom_within_bounds + num_of_chrom_beyond
            )
            length_stats = ChromLengthStats(
                oobr=sensitivity,
                beyond_range=chroms_beyond_range,
                num_of_chrom_beyond=num_of_chrom_beyond,
                percentage_bed_chrom_beyond=(
                    100 * num_of_chrom_beyond / len(bed_chrom_set)
                ),
                percentage_genome_chrom_beyond=(
                    100 * num_of_chrom_beyond / len(genome_chrom_set)
                ),
            )

        else:
            length_stats = ChromLengthStats()

        # Layer 3 Calculate Sequence Fit if any query chrom names were present
        if len(query_keys_present) > 0:
            bed_sum = 0
            ref_genome_sum = 0
            for q_chr in query_keys_present:
                bed_sum += int(genome_chrom_sizes[q_chr])
            for g_chr in genome_chrom_sizes:
                ref_genome_sum += int(genome_chrom_sizes[g_chr])

            seq_fit_stats = SequenceFitStats(sequence_fit=bed_sum / ref_genome_sum)

        else:
            seq_fit_stats = SequenceFitStats(sequence_fit=None)

        return CompatibilityStats(
            chrom_name_stats=name_stats,
            chrom_length_stats=length_stats,
            chrom_sequence_fit_stats=seq_fit_stats,
        )

    def get_igd_overlaps(self, bedfile: str) -> Union[dict[str, dict], dict[str, None]]:
        """
        This is the third layer compatibility check.
        It runs helper functions to execute an igd search query across an Integrated Genome Database.

        It returns a dict of dicts containing keys (file names) and values (number of overlaps). Or if no overlaps are found,
        it returns an empty dict.

        Currently for this function to work, the user must install the C version of IGD and have created a local igd file
        for the Excluded Ranges Bedset:
        https://github.com/databio/IGD
        https://bedbase.org/bedset/excluderanges

        """
        if not self.igd_path:
            return {"igd_stats": None}
        # Construct an IGD command to run as subprocess
        igd_command = IGD_LOCATION + f" search {self.igd_path} -q {bedfile}"

        returned_stdout = run_igd_command(igd_command)

        # print(f"DEBUG: {returned_stdout}")

        if not returned_stdout:
            return {"igd_stats": None}

        igd_overlap_data = parse_IGD_output(returned_stdout)

        if not igd_overlap_data:
            return {
                "igd_stats": {}
            }  # None tells us if the bed file never made it to layer 3 or perhaps igd errord, empty dict tells us that there were no overlaps found
        else:
            overlaps_dict = {}
            for datum in igd_overlap_data:
                if "file_name" in datum and "number_of_hits" in datum:
                    overlaps_dict.update({datum["file_name"]: datum["number_of_hits"]})

        return overlaps_dict

    def determine_compatibility(
        self,
        bedfile: str,
        ref_filter: Optional[List[str]] = None,
        concise: Optional[bool] = False,
    ) -> Union[Dict[str, CompatibilityStats], Dict[str, CompatibilityConcise]]:
        """
        Given a bedfile, determine compatibility with reference genomes (GenomeModels) created at Validator initialization.

        :param str bedfile: path to bedfile on disk, .bed
        :param list[str] ref_filter: list of ref genome aliases to filter on.
        :param bool concise: if True, only return a concise list of compatibility stats
        :return list[dict]: a list of dictionaries where each element of the array represents a compatibility dictionary
                            for each refgenome model.
        """

        if ref_filter:
            # Before proceeding filter out unwanted reference genomes to assess
            for genome_model in self.genome_models:
                if genome_model.genome_alias in ref_filter:
                    self.genome_models.remove(genome_model)

        bed_chrom_info = get_bed_chrom_info(
            bedfile
        )  # for this bed file determine the chromosome lengths

        if not bed_chrom_info:
            # if there is trouble parsing the bed file, return None
            raise ValidatorException

        model_compat_stats = {}

        for genome_model in self.genome_models:
            # First and Second Layer of Compatibility
            model_compat_stats[genome_model.genome_alias]: CompatibilityStats = (
                self.calculate_chrom_stats(bed_chrom_info, genome_model.chrom_sizes)
            )
            # Fourth Layer is to run IGD but only if layer 1 and layer 2 have passed
            if (
                model_compat_stats[
                    genome_model.genome_alias
                ].chrom_name_stats.passed_chrom_names
                and not model_compat_stats[
                    genome_model.genome_alias
                ].chrom_length_stats.beyond_range
            ):
                model_compat_stats[genome_model.genome_alias].igd_stats = (
                    self.get_igd_overlaps(bedfile)
                )

            # Once all stats are collected, process them and add compatibility rating
            model_compat_stats[genome_model.genome_alias].compatibility = (
                self.calculate_rating(model_compat_stats[genome_model.genome_alias])
            )

        if concise:
            concise_dict = {}
            for name, stats in model_compat_stats.items():
                concise_dict[name] = self._create_concise_output(stats)

            return concise_dict

        return model_compat_stats

    def calculate_rating(self, compat_stats: CompatibilityStats) -> RatingModel:
        """
        Given compatibility stats for a specific ref genome determine the compatibility tier

        Tiers Definition ----

        Tier1: Excellent compatibility, 0 pts
        Tier2: Good compatibility, may need some processing, 1-3 pts
        Tier3: Bed file needs processing to work (shifted hg38 to hg19?), 4-6 pts
        Tier4: Poor compatibility, 7-9 pts

        :param CompatibilityStats compat_stats: dicitionary containing unprocessed compat stats
        :return dict : containing both the original stats and the final compatibility rating for the reference genome
        """

        points_rating = 0

        # 1. Check extra sequences sensitivity and assign points based on how sensitive the outcome is
        # sensitivity = 1 is considered great and no points should be assigned
        xs = compat_stats.chrom_name_stats.xs
        if xs < 0.3:
            points_rating += 6  # 3 + 1 + 1 + 1
        elif xs < 0.5:
            points_rating += 5  # 3 + 1 + 1
        elif xs < 0.7:
            points_rating += 4  # 3 + 1
        elif xs < 1:
            points_rating += 3
        else:
            pass

        # 2. Check OOBR and assign points based on sensitivity
        # only assessed if no extra chroms in query bed file
        if compat_stats.chrom_name_stats.passed_chrom_names:
            oobr = compat_stats.chrom_length_stats.oobr

            if oobr < 0.3:
                points_rating += 6  # 3 + 1 + 1 + 1
            elif oobr < 0.5:
                points_rating += 5  # 3 + 1 + 1
            elif oobr < 0.7:
                points_rating += 4  # 3 + 1
            elif oobr < 1:
                points_rating += 3
        else:
            # Do nothing here, points have already been added when Assessing XS if it is not == 1
            pass

        # Check Sequence Fit - comparing lengths in queries vs lengths of queries in ref genome vs not in ref genome
        sequence_fit = compat_stats.chrom_sequence_fit_stats.sequence_fit
        if sequence_fit:
            # since this is only on keys present in both, ratio should always be less than 1
            # Should files be penalized here or actually awarded but only if the fit is really good?
            if sequence_fit < 0.90:
                points_rating += 1
            if sequence_fit < 0.60:
                points_rating += 1
            if sequence_fit < 0.60:
                points_rating += 1

        else:
            # if no chrom names were found during assessment
            points_rating += 4

        # Run analysis on igd_stats
        # WIP, currently only showing IGD stats for informational purposes
        if compat_stats.igd_stats and compat_stats.igd_stats != {}:
            self.process_igd_stats(compat_stats.igd_stats)

        tier_ranking = 0
        if points_rating == 0:
            tier_ranking = 1
        elif 1 <= points_rating <= 3:
            tier_ranking = 2
        elif 4 <= points_rating <= 6:
            tier_ranking = 3
        elif 7 <= points_rating:
            tier_ranking = 4
        else:
            print(
                f"Catching points discrepancy,points = {points_rating}, assigning to Tier 4"
            )
            tier_ranking = 4

        return RatingModel(assigned_points=points_rating, tier_ranking=tier_ranking)

    def process_igd_stats(self, igd_stats: dict):
        """
        Placeholder to process IGD Stats and determine if it should impact tier rating
        """
        pass

    @staticmethod
    def _build_default_models() -> list[GenomeModel]:
        """
        Builds a default list of GenomeModels from the chrom.sizes folder.
        Uses file names as genome alias.

        return list[GenomeModel]
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))

        chrm_sizes_directory = os.path.join(dir_path, "chrom_sizes")

        all_genome_models = []
        for root, dirs, files in os.walk(chrm_sizes_directory):
            for file in files:
                if file.endswith(".sizes"):
                    curr_genome_model = GenomeModel(
                        genome_alias=file, chrom_sizes_file=os.path.join(root, file)
                    )
                    all_genome_models.append(curr_genome_model)

        return all_genome_models

    @staticmethod
    def _create_concise_output(output: CompatibilityStats) -> CompatibilityConcise:
        return CompatibilityConcise(
            xs=output.chrom_name_stats.xs,
            oobr=output.chrom_length_stats.oobr,
            sequence_fit=output.chrom_sequence_fit_stats.sequence_fit,
            assigned_points=output.compatibility.assigned_points,
            tier_ranking=output.compatibility.tier_ranking,
        )


def run_igd_command(command):
    """Run IGD via a subprocess, this is a temp implementation until Rust IGD python bindings are finished."""

    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        return result.stdout
    else:
        print(f"Error running command: {result.stderr}")
        return None


def parse_IGD_output(output_str) -> Union[None, List[dict]]:
    """
    Parses IGD terminal output into a list of dicts
    Args:
      output_str: The output string from IGD

    Returns:
      A list of dictionaries, where each dictionary represents a record.
    """

    try:
        lines = output_str.splitlines()
        data = []
        for line in lines:
            if line.startswith("index"):
                continue  # Skip the header line
            elif line.startswith("Total"):
                break  # Stop parsing after the "Total" line
            else:
                fields = line.split()
                record = {
                    "index": int(fields[0]),
                    "number_of_regions": int(fields[1]),
                    "number_of_hits": int(fields[2]),
                    "file_name": fields[3],
                }
                data.append(record)
        return data
    except Exception:
        return None
