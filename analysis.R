# Analysis
library(tidyverse)

uniq <- function(x) {
  x %>% 
    nest(d = -uid) %>% 
    mutate(d = map(
      d, 
      function(d) {
        if (nrow(d) == 1) return(d)
        d %>% 
          mutate(err = ws_error + tu_error) %>% 
          nest(a = -err) %>%
          filter(err == min(err)) %>%
          unnest(cols = a) %>%
          select(-err) %>% 
          arrange(abs(ws_end) + abs(tu_end)) %>%
          .[1,]
      }
    )) %>%
    unnest(cols = d)
}

## Load the old data
load('data/thesis-parameter-estimation.rda')
R <- recovered_parameters %>% mutate(method = 'old') %>% uniq()
S <- status %>% mutate(method = 'old')

## Load the new data with untethered step sizes
load('data/thesis-parameter-estimation-uncoupled.rda')
R <- bind_rows(
  R, 
  mutate(recovered_parameters, method = 'uncoupled') %>% uniq()
)
S <- bind_rows(S, mutate(status, method = 'uncoupled'))

## Load the new data with untethered step sizes and fixed TU
load('data/thesis-parameter-estimation-uncoupled-fixedTU.rda')
R <- bind_rows(
  R, 
  mutate(recovered_parameters, method = 'fixedTU') %>% uniq()
)
S <- bind_rows(S, mutate(status, method = 'fixedTU'))

# Check timings
# Remember the old data include shuffles, so will be ~10x longer anyway
S %>% 
  mutate(duration = time_end - time_start) %>%
  ggplot(aes(duration * cores, fill = method)) +
  geom_histogram(bins = 100, position = 'identity', alpha = .5)

# Check differences in parameters
R_diff <- R %>%
  select(uid, method, ends_with('_error'), ends_with('_end')) %>%
  pivot_longer(c(-uid, -method)) %>%
  pivot_wider(names_from = method) %>%
  mutate(diff_uncoupled = old - new, diff_fixedTU = old - fixedTU) %>%
  select(c(-old, -new, -fixedTU))

R_diff %>%
  pivot_longer(
    starts_with('diff_'), 
    names_to = "diff", 
    names_pattern = "diff_(.+)"
  ) %>%
  ggplot(aes(value, fill = diff)) +
  geom_histogram(bins = 100, alpha = 1/3, position = 'identity') +
  facet_wrap(~name, scales = 'free')
# So most participants' values are stable across implementations, although some
# have some quite dramatic differences. Notably, some participants have much
# higher new values for weighted selection than others, though plenty of others
# have lower ones.

for (v in c('ws_end', 'tu_end', 'ws_error', 'tu_error')) {
  print(
    R %>%
      select(method, !!v) %>%
      pivot_longer(-method) %>%
      ggplot(aes(x = value)) +
      geom_histogram(bins = 100) +
      facet_grid(method~name, scales = 'free_x', as.table = F)
  )
}


R %>%
  group_by(method) %>%
  select(method, ends_with('_end')) %>%
  summarise(across(.fns = list(median = median, mean = mean, sd = sd))) %>% 
  mutate(across(.cols = -method, ~ round(., 3)))

R %>% 
  group_by(method) %>%
  select(method, ws_end) %>%
  nest(d = -method) %>%
  mutate(
    d = map(
      d, 
      ~ PearsonDS::empMoments(.$ws_end) %>%
        as.data.frame() %>%
        t() %>%
        as_tibble()
    )
  ) %>%
  unnest(cols = d) %>%
  transmute(method, mean, sd = sqrt(variance), skewness, kurtosis)
