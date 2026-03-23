## How to reconstruct signal without phase shift?

set.seed(1234)
simdat <- tibble(time = 0:(499),
                 c1 = 1.5 * sin(2 * pi * 1/(12) * time + 0),
                 c2 = .5 * sin(2 * pi * 1/(10) * time + pi/2),
                 c3 = 1 * sin(2 * pi * 1/30 * time + pi/3),
                 y = c1 + c2 + c3,
                 r = rnorm(500, 0, 1),
                 x = y + r)

test_comp <- matrix(nrow = nrow(simdat), ncol = 4)
test_f <- c(0, 1/12, 1/10, 1/30)

for (i in 1:ncol(test_comp)) {
  test_comp[, i] <- kzft(simdat$x, f = test_f[i], m = 120, k = 2)
}

## reconstruct by 2 * sum of kzft's and extracting real part
recon <- tibble(t = simdat$time,
                x = simdat$x,
                y = simdat$y,
                cplx_signal = rowSums(test_comp),
                re_signal = 2*Re(cplx_signal))

recon %>%
  filter(t < 200) %>%
  ggplot(aes()) +
  geom_line(aes(t, y)) +
  # geom_line(aes(t, re_signal), color = "red") +
  geom_line(aes(t+2, re_signal), color = "blue")



#### 
n <- nrow(test_comp)
c <- ncol(test_comp)
recon <- matrix(nrow = n, ncol = c)

for (i in 1:c) {
  for (j in 1:n) {
    recon[j,i] <- 1/n * sum(test_comp[,i] * exp(complex(real = 0, imaginary = 2*pi/n*(j-1) * (0:(n-1)))))
  }
}

head(recon)



tibble(t = simdat$time, x = simdat$x,
       cplx = rowSums(test_comp)) %>%
  mutate(xn = sum(cplx * exp(complex(real = 0, imaginary = 2 * pi * t))))
       
       y = simdat$y) %>% 
  filter(t < 200)

tibble(t = simdat$time,
       x = simdat$x,
       y = simdat$y,
       pred = 2*Re(rowSums(test_comp)),
       theta = Arg(rowSums(test_comp)),
       p = pred * cos()) %>% head
  filter(t < 200) %>%
  ggplot(aes()) +
  geom_line(aes(t, x)) +
  geom_line(aes(t, y), color = "red") +
  geom_line(aes(t, pred), color = "blue")



