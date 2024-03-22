test_that("moving_average works", {
    expect_equal(
        moving_average(1:5, 3),
        c(2, 2, 2, 3, 4)
    )
})
