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
  
  
  
  ## eda
  dat %>%
    mutate(kz1 = kza::kz(.$x, m = 24*7*3+1, k = 3),
           kz2 = kza::kz(.$x, m = 24*7*5+1, k = 3),
           kz3 = kza::kz(.$x, m = 24*7*7+1, k = 3),
           kz4 = kza::kz(.$x, m = 24*7*11+1, k = 3)) %>%
    select(t_start, kz1, kz2, kz3, kz4) %>%
    tidyr::pivot_longer(cols = c(kz1, kz2, kz3, kz4),
                        values_to = "energy", names_to = "name") %>%
    ggplot(aes(t_start, energy, color = name)) +
    geom_line()
  
  ## estimate key constituents over years
  dat2 <- tibble(t_start = seq(from = as.POSIXct("2016-03-24 00:00:00", tz = 'EST'),
                               to = as.POSIXct("2023-04-01 00:00:00", tz = 'EST'),
                               by = "hour"),
                 t_end = t_start + 60^2) %>%
    mutate(t = row_number()) %>%
    relocate(t) %>%
    left_join(cap, join_by(overlaps(t_start, t_end, Start_Time, End_Time, bounds = "[)") )) %>%
    mutate(energy = difftime(pmin(t_end, End_Time),
                             pmax(t_start, Start_Time), units = "hours") %>%
             as.numeric() %>%
             tidyr::replace_na(0)) %>%
    filter(t_start >= as.POSIXct("2019-01-01 00:00:00", tz = 'EST'),
           t_start < as.POSIXct("2023-03-01 00:00:00", tz = 'EST')) %>%
    group_by(t, t_start, t_end) %>%
    summarise(x = sum(energy)) %>%
    ungroup()
  
  f_components2 <- matrix(nrow = nrow(dat2), ncol = nrow(f_include))
  
  for (i in 1:ncol(f_components2)) {
    f_components2[, i] <- kzft(dat2$x, f = f_include$f[i], m = 168*13, k = 5)
  }
  
  
  dat2 %>%
    filter(t_start >= as.POSIXct("2021-01-01 00:00:00", tz = 'EST'),
           t_start < as.POSIXct("2022-01-01 00:00:00", tz = 'EST')) %>%
    mutate(Day = lubridate::wday(t_start, label = TRUE) %>%
             factor(levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")),
           Hour = hour(t_start)) %>%
    group_by(Day, Hour) %>%
    summarise(x = sum(x)) %>%
    ungroup() %>%
    mutate(rel_risk = x / mean(x)) %>%
    ggplot(aes(x = Hour, y = rel_risk)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 1) +
    theme(panel.spacing.x = unit(0, "lines"),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = c(0, 8, 17),
                       expand = c(0,0),
                       name = "Day of Week - Hour") +
    scale_y_continuous(name = "Total Energy",
                       labels = scales::label_number(accuracy = .1)) +
    facet_grid(~ Day)


