test_that(
    "bucket_class works",
    {
        expect_equal(
            calc_cn_runoff(list(CN = 1.5)),
            1.5
        )
    }
)