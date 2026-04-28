library(tidyverse)

ggplot(BOD, aes(Time, demand))+
  geom_point(seize=3)+
  geom_line(color="forestgreen")

CO2 %>% ggplot(aes(conc, uptake, colour = Treatment)) +
  geom_point(size=3, alpha=0.5) +
  geom_smooth(se = F)+
  facet_wrap(~Type)+
  labs(title = 'Concentration of CO2') + 
  theme_bw()

CO2 %>% ggplot(aes(Treatment, uptake))+    # %>% operator (pipe operator) directly gives in input CO2 to the data= parameter
  geom_boxplot()+   # geometric display
  geom_point(alpha=0.5,
    aes(size=conc, color=Plant))+ # size by concentration, color by different plant
  coord_flip()+  # flip the coordinates 
  theme_bw()+   # many to choose from
  labs(title='Uptake upon treatment')+
  facet_wrap(~Type)  # we split the plot in two columns, according to the formula 

mpg %>%
  filter(cty<25)%>%
  ggplot(aes(displ, cty))+
  geom_point(aes(colour = drv, size = trans),
             alpha=0.5)+
  geom_smooth(method = 'glm')+
  facet_wrap(~year, nrow=1)+
  labs(x='Engine size',
       y='MPG in the city',
       title = 'Fuel efficiency')+
  theme_bw()

msleep %>%
  drop_na(vore)%>%
  ggplot(aes(fct_infreq(vore)))+
  geom_bar(fill = 'darkcyan')+
  theme_bw()+
  labs(x='Who eats what?',
       y='Frequency',
       title = 'Number of observation per order')

msleep %>%
  filter(bodywt < 5) %>%
  ggplot(aes(bodywt, brainwt))+
  geom_point(aes(colour = sleep_total, size=awake))+
  geom_smooth(method = lm)+
  labs(x='Body weight',
       y='Brain weight',
       title = 'Brain vs body weight')+
  theme_bw()

msleep %>% 
  drop_na(vore, sleep_rem) %>%
  group_by(vore) %>% 
  summarise('Average total sleep' = mean(sleep_total),
            'Max. REM sleep' = max(sleep_rem))
