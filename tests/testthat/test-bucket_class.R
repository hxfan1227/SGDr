# Test the Bucket class
test_that("Bucket class getters return correct values", {
  # Test case 1
  bucketParams <- list(
    layer = 1,
    WP = 100,
    z = 200,
    FCmm = 300,
    SAT = 400,
    Swres = 500,
    Ksat = 600,
    n = 0.5,
    phi_soil = 0.3,
    Y = 0.2,
    Z = 0.1
  )
  bucket <- new(Bucket, bucketParams)
  expect_equal(bucket$get_WP(), 100)
  expect_equal(bucket$get_z(), 200)
  expect_equal(bucket$get_FCmm(), 300)
  expect_equal(bucket$get_SAT(), 400)
  expect_equal(bucket$get_Swres(), 500)
  expect_equal(bucket$get_Ksat(), 600)
  expect_equal(bucket$get_n(), 0.5)
  expect_equal(bucket$get_phi_soil(), 0.3)
  expect_equal(bucket$get_Y(), 0.2)
  expect_equal(bucket$get_Z(), 0.1)

  # Test case 2
    bucketParams <- list(
        layer = 2,
        WP = 200,
        z = 300,
        FCmm = 400,
        SAT = 500,
        Swres = 600,
        Ksat = 700,
        n = 0.6,
        phi_soil = 0.4
    )
    bucket <- new(Bucket, bucketParams)
    expect_equal(bucket$get_WP(), 200)
    expect_equal(bucket$get_z(), 300)
    expect_equal(bucket$get_FCmm(), 400)
    expect_equal(bucket$get_SAT(), 500)
    expect_equal(bucket$get_Swres(), 600)
    expect_equal(bucket$get_Ksat(), 700)
    expect_equal(bucket$get_n(), 0.6)
    expect_equal(bucket$get_phi_soil(), 0.4)
    expect_equal(bucket$get_Y(), 0)
    expect_equal(bucket$get_Z(), 0)
# Test case 3
    bucketParams <- list(
        layer = 1,
        WP = 300,
        z = 400,
        FCmm = 500,
        SAT = 600,
        Swres = 700,
        Ksat = 800,
        n = 0.7,
        phi_soil = 0.5
    )
    expect_error(new(Bucket, bucketParams))
})