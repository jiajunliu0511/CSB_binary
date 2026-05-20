library(tidyverse)
library(scales)
library(cowplot)
library(latex2exp)
theme_set(theme_bw())
set.seed(2024)


# est ---------------------------------------------------------

res_est <- map_dfr(1:4, function(case) {
  load(str_glue("output/cache/inf_{case}.RData"))
  load("data/setup.RData")
  inf_result %>% 
    left_join(setup, by ="case") %>% 
    filter(!str_detect(method, "FRT")) %>% 
    mutate(
      SM = ifelse(sm, "Correct", "Wrong"),
      OM = ifelse(om, "Correct", "Wrong"),
      Method = case_when(
        method == "No Borrow DiM" ~ "No Borrow DifMean",
        method == "No Borrow AIPW" ~ "No Borrow CovAdj",
        TRUE ~ method
      ) %>% 
        as_factor() %>% 
        fct_relevel("No Borrow DifMean","No Borrow CovAdj",
                    "Borrow Naive", "Borrow IPW", "Borrow CW"),
      `tauhat - tau` = est - mean(tau)
    ) %>% 
    arrange(rep, Method) %>% 
    select(SM, OM, Method, `tauhat - tau`)
}) %>% mutate(across(c("SM", "OM"), as_factor))

mse_labels <- res_est %>%
  group_by(SM, OM, Method) %>%
  summarise(MSE = mean((`tauhat - tau`)^2, na.rm = T)) %>%
  mutate(
    RelMSE = round(MSE / MSE[1] * 100),
    rank_MSE = rank(MSE),
    is_min = rank_MSE <= 2,
    label_type = "Relative MSE"
  ) %>% 
  ungroup()

res_est %>% 
  ggplot(aes(Method, `tauhat - tau`, color = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~SM + OM, nrow = 1,labeller = label_both) +
  ylab(TeX("\\hat{\\tau} - \\tau")) +
  theme(axis.text.x = element_blank()) +
  geom_label(
    data = mse_labels,
    aes(x = Method, y = max(res_est$`tauhat - tau`, na.rm = T) + 0.1, 
        label = RelMSE,
        fontface = ifelse(is_min, "bold", "plain"),
        fill = label_type),
    inherit.aes = FALSE,
    color = "black",  
    vjust = 0.8,
    size = 2.8
  ) +
  scale_fill_manual(values = c("Relative MSE" = "white")) +
  labs(fill = "", color = "") +
  theme(legend.position="bottom",legend.box = "vertical") +
  guides(color = guide_legend(nrow = 1, order = 2),
         fill = guide_legend(nrow = 1, order = 1))

ggsave("chart/sim_est.pdf", width = 10, height = 6)

# type I ------------------------------------------------------------------

res_test <- map_dfr(1:4, function(case) {
  load(str_glue("output/metrics_{case}.RData"))
  load("data/setup.RData")
  metrics %>% 
    left_join(setup, by ="case") %>% 
    mutate(
      SM = ifelse(sm, "Correct", "Wrong"),
      OM = ifelse(om, "Correct", "Wrong"),
      Test = ifelse(
        str_detect(method, "FRT"),
        "FRT",
        "Asymptotic"
      ) %>% as_factor,
      method = gsub("\\+FRT", "", method),
      Method = case_when(
        method == "No Borrow DiM" ~ "No Borrow DifMean",
        method == "No Borrow AIPW" ~ "No Borrow CovAdj",
        TRUE ~ method
      ) %>% 
        as_factor() %>% 
        fct_relevel("No Borrow DifMean", "No Borrow CovAdj",
                    "Borrow Naive", "Borrow IPW", "Borrow CW")
    ) %>% 
    arrange(Method, Test) %>% 
    select(SM, OM, Method, Test, `Type I`, Power)
}) %>% mutate(across(c("SM", "OM"), as_factor))


p1 <- res_test %>% 
  ggplot(aes(Method, `Type I`, color = Method)) +
  geom_hline(yintercept = 0.05, linetype = 1) +
  geom_point(aes(shape = Test), size = 2, stroke = 1) +
  facet_wrap(SM ~ OM, nrow = 1, labeller = label_both) +
  ylab(TeX("Type I Error Rate")) +
  theme(axis.text.x = element_blank()) +
  #ylim(0, 0.2) +
  scale_shape_manual(values = c(3, 16)) +
  labs(shape= "", color = "") +
  theme(legend.position="bottom",legend.box = "vertical") +
  guides(color = guide_legend(nrow = 1, order = 2),
         shape = guide_legend(nrow = 1, order = 1))


p2 <- res_test %>% 
  ggplot(aes(Method, Power, color = Method)) +
  geom_point(aes(shape = Test), size = 2, stroke = 1) +
  facet_wrap(SM ~ OM, nrow = 1, labeller = label_both) +
  theme(axis.text.x = element_blank()) +
  scale_shape_manual(values = c(3, 16)) +
  labs(shape= "", color = "") +
  theme(legend.position="bottom",legend.box = "vertical") +
  guides(color = guide_legend(nrow = 1, order = 2),
         shape = guide_legend(nrow = 1, order = 1))


legend <- get_plot_component(p1, "guide-box", return_all = T)[[3]]

combined_plot <- plot_grid(
  p1+ theme(legend.position="none"), 
  p2+ theme(legend.position="none"),
  nrow = 2,
  align = "v"
)


plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.2))


ggsave("chart/sim_test.pdf", width = 10, height = 6)


