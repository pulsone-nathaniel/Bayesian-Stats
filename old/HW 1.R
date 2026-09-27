# Quesiton 1

y = read.delim("algae.txt")
y = 44
n = 274

pbeta(.2, 46, 240)
qbeta(.95, 46, 240)

theta0 = seq(0.1, 0.9, 0.1)
n0 = c(1, 2, 8, 16, 32)

Epost = matrix(nrow = 9, ncol = 5)
for(i in 1:9){
  for(j in 1:5){
    w = n0[j]
    Epost[i, j] = (1 / (n + w))*y + (w / (n + w))*theta0[i]
  }
}
contour(theta0, n0, Epost, main = "Posterior Expectation",
        ylab = "Prior Sample Size", xlab = "Prior Theta")

# Question 2
sa = sum(c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6))
na = 10

sb = sum(c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7))
nb = 13
qgamma(.025, 125, 14)

n0 = seq(1, 50, 1)
a = 12*n0 + sb
b = n0 + 12
a/b
b

# Question 4
q = seq(0, 1, 0.01)
pt = density(rbeta(200000, 8, 2))
pyt = 151532656696 * q^15 * (1-q)^28
pty = density(rbeta(200000, 23, 30))
plot(pty, type = "l", main = "Posterior Density")
SSs
