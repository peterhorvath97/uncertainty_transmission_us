estout_clean <- read_rds('US_only/data/estout_clean.Rds')

irfs <- estout_clean[['irfs']]


estout_clean$invalids %>% 
  mutate(index = case_when(index == 'epu_debt' ~ 'EPU: Debt, crises',
                           index == 'epu_finreg' ~ 'EPU: Fin. Regulation',
                           index == 'epu_fis' ~ 'EPU: Fiscal Policy',
                           index == 'epu_mon' ~ 'EPU: Monetary Policy',
                           index == 'epu_reg' ~ 'EPU: Regulation',
                           index == 'epu_trade' ~ 'EPU: Trade Policy'),
         threshold = case_when(threshold == '0.5' ~ 'Median',
                               threshold == '0.7' ~ '70th Perc.'),
         control = case_when(control == 'none' ~ 'No controls',
                             control == 'uncertainty' ~ 'Uncertainty',
                             control == 'financial_stress' ~ 'Financial Stress',
                             control == 'both' ~ 'Uncertainty and Financial Stress'),
         spec = case_when(spec == 'main' ~ 'Interest rate measured by FFR',
                          spec == 'lag' ~ 'Regime indicator from lagged level of uncertainty',
                          spec == 'r1y'  ~ 'Interest rate measured by')) %>% 
  janitor::tabyl(spec)

lw = .75




  


irfs %>% 
  select(-mult, -restr) %>% 
  filter(label != 'R',
         index != 'epu') %>% 
  filter(spec == 'shadow',
         control == 'none',
         thresh == '0.5') %>% 
  mutate(index = case_when(index == 'epu_debt' ~ 'EPU: Debt, crises',
                           index == 'epu_finreg' ~ 'EPU: Fin. Regulation',
                           index == 'epu_fis' ~ 'EPU: Fiscal Policy',
                           index == 'epu_mon' ~ 'EPU: Monetary Policy',
                           index == 'epu_reg' ~ 'EPU: Regulation',
                           index == 'epu_trade' ~ 'EPU: Trade Policy')) %>% 
  ggplot(aes(x = t, y = median, color = regime
  )) +
  geom_line(linewidth = lw) +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(label~index, scales = "free", nrow = 2)+
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y ="")





  