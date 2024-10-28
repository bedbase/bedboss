def runn():
    import cProfile
    import pstats

    from bedboss.bedboss import run_all

    # with cProfile.Profile() as pr:
    run_all(
        bedbase_config="/home/bnt4me/virginia/repos/bbuploader/config_db_local.yaml",
        outfolder="/home/bnt4me/virginia/repos/bbuploader/data",
        genome="hg38",
        # input_file="/home/bnt4me/virginia/repos/bedboss/test/data/bed/hg38/GSM6732293_Con_liver-IP2.bed",
        input_file="/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/49666564c051c7ece0cacff6b74e1d90.bed.gz",
        # input_file="/home/bnt4me/virginia/repos/bedboss/scripts/stats/data/b72912f88165f2f252f273f4c3043fe3.bed.gz",
        # input_file="/home/bnt4me/virginia/repos/bedboss/scripts/stats/8f37f6b696eb1862328efe7237ddccdf.bed.gz",
        input_type="bed",
        force_overwrite=True,
        upload_pephub=True,
        upload_s3=True,
        upload_qdrant=True,
        name="test",
    )

    # stats = pstats.Stats(pr)
    # stats.sort_stats(pstats.SortKey.TIME)
    # stats.dump_stats(filename="test_profiling")


if __name__ == "__main__":
    runn()
