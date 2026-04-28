library(tidyverse)
library(naniar)
library(gtExtras)

# Magnitude of missing data ####

## Table ####
miss_var_summary(airquality) %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Missingness of variables')

## Plot ####
gg_miss_var(airquality)


# Distribution of missing data ####

## Table of observations with missing data ####
airquality %>% 
  filter(!complete.cases(.)) %>%
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Rows that contain missing data')

## Heatmap of missing data ####
vis_miss(airquality)


# Relationship to one variable ####

airquality %>% 
  mutate(missing_ozone = factor(is.na(Ozone),
                                levels = c('TRUE', 'FALSE'),
                                labels = c('Missing', 'Not Missing'))) %>% 
  ggplot(aes(x = Wind, fill = missing_ozone))+
  geom_histogram(position = 'stack')+
  labs(title = 'Distribution of Wind speeds for Missing vs. Non-Missing',
       x = 'Wind speed',
       y = 'Ozone observation',
       fill = 'Missingness')+
  theme_bw()

airquality %>% 
  mutate(missing_ozone = factor(is.na(Ozone),
                                levels = c('FALSE', 'TRUE'),
                                labels = c('Not Missing', 'Missing'))) %>% 
  ggplot(aes(x = Wind, fill = missing_ozone))+
  geom_histogram(position = 'stack')+
  labs(title = 'Distribution of Wind speeds for Missing vs. Non-Missing',
       x = 'Wind speed',
       y = 'Ozone observation',
       fill = 'Missingness')+
  theme_bw()

## or more variables ####
airquality %>% 
  select(Ozone, Solar.R, Wind, Temp) %>% 
  ggplot(aes(x = Wind,
             y = Temp,
             size = Solar.R,
             colour = is.na(Ozone)))+
  geom_point(alpha=0.6)+
  facet_wrap(~is.na(Ozone))+
  labs(title = 'Missing Ozone data by Wind and Temperature',
       x = 'Wind speed',
       y = 'Temperature',
       color = 'Missing Ozone data')+
  theme_bw()


# Dealing with missing values ####

## Drop it ####
starwars %>% 
  select(name, mass, height, hair_color) %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Starwars charachters') %>% 
  gt_highlight_rows(rows = is.na(mass), fill = 'steelblue') %>% 
  gt_highlight_rows(rows = is.na(hair_color), fill = 'lightpink')
starwars %>% 
  select(name, mass, height, hair_color) %>% 
  drop_na(mass) %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Starwars charachters')

## Change it ####
starwars %>% 
  select(name, hair_color, species) %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Starwars charachters')
starwars %>% 
  select(name, hair_color, species) %>%
  filter(species == 'Droid') %>% 
  mutate(hair_color = case_when(is.na(hair_color) & species == 'Droid'
                                ~ 'none', TRUE ~ hair_color)) %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Starwars charachters')

## Impute it ####
starwars %>% 
  mutate(height = case_when(
    is.na(height) ~ median(starwars$height, na.rm = TRUE),  # if height has a NA put the median of the dataset
  TRUE ~ height)) %>%   # if previous conditions are not met, leave the value as it was
  select(name, height) %>% 
  arrange(name) %>% 
  gt() %>% 
  gt_theme_nytimes() %>% 
  tab_header(title = 'Starwars charachters') %>% 
  gt_highlight_rows(rows = name %in% c('Arvel Crynyd', 'BB8', 'Poe Dameron', 'Captain Phasma'))
