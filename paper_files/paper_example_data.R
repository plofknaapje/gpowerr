
data <-
  data.frame(
    a = 1:6,
    b = 1:6 ^ 2,
    c = c(0,-1, 0.5, 0,-2,-3),
    d = c(1, 1, 2, 2, 3, 3),
    e = c(0, 1, 0, 1, 0, 1)
  )
gpower(data, 2, 0.1)
prcomp(data)
