library(tidyverse)

starwars %>%
  select(gender, mass, height, species) %>%
  filter(species=='Human') %>%
  na.omit()%>%
  mutate(height=height/100)%>%
  mutate(BMI=mass / height^2)%>%
  group_by(gender)%>%
  summarise(Average_BMI=mean(BMI))

# rename variables
starwars %>%
  rename(home = homeworld)%>%
  glimpse()

# reorder variables
starwars %>%
  select(sex, gender, everything())%>%
  glimpse()

# change variable type
starwars %>%
  mutate(gender=as.factor(gender)) %>%
  glimpse()

# select variables to work with
starwars %>%
  select(2:4,
         species,
         starts_with('veh'),
         contains('ship'))%>%
  glimpse()

# filter and arrange data
starwars %>%
  filter((species == 'Human' | species == 'Droid') & height > 165) %>%
  select(species, height, gender, hair_color) %>%
  arrange(height) %>%
  glimpse()
# or
starwars %>%
  filter(species %in% c('Human', 'Droid') & height > 165) %>%
  select(species, height, gender, hair_color) %>%
  arrange(height) %>%
  glimpse()

# conditional changes
starwars %>%
  select(species, height, starships, name) %>%
  drop_na(starships) %>%
  mutate(size = if_else(height > 180, 'big', 'small', missing = 'unknown')) %>%
  glimpse()

# recode data and rename a variable
starwars %>%
  select(species, height, starships, name, gender) %>%
  drop_na(starships) %>%
  mutate(size = if_else(height > 180, 'big', 'small', missing = 'unknown')) %>%
  mutate(size = recode(size, 'big' = 1, 'small' = 0)) %>%
  mutate(gender = recode(gender, male = 'm',
                         female = 'f')) %>% 
  glimpse()
  
# reshape dataset from long to wide (or viceversa)
library(gapminder)
view(gapminder)
data <- select(gapminder, country, year, lifeExp)
view(data)   # now it is long
wide_data <- data %>%
  pivot_wider(names_from = year, values_from = lifeExp)
view(wide_data)  # now it is wide

long_data <- wide_data %>%
  pivot_longer(cols=2:13,
               names_to = 'year',
               values_to = 'lifeExp')
view(long_data)

# summarise
msleep %>% 
  drop_na(vore) %>% 
  group_by(vore) %>% 
  summarise(Lower = min(sleep_total),
            Average = mean(sleep_total),
            Upper = max(sleep_total),
            Difference = max(sleep_total)-min(sleep_total)) %>% 
  arrange(Average) %>% 
  view()


# Factors ####
library(forcats)
View(gss_cat)

## Manually reordering ####
gss_cat %>% 
  pull(race) %>% 
  levels()

gss_cat %>% 
  select(race) %>% 
  table()

gss_cat %>% 
  mutate(race = fct_drop(race)) %>% 
  pull(race) %>% 
  levels()

gss_cat %>% 
  mutate(race = fct_drop(race)) %>% 
  mutate(race = fct_relevel(race,
                            c("White", "Black", "Other"))) %>% 
  pull(race) %>% 
  levels()

## By the average value of another variable ####
gss_cat %>% 
  drop_na(tvhours) %>% 
  group_by(relig) %>% 
  summarise(mean_tv = mean(tvhours)) %>% 
  ggplot(aes(mean_tv, relig)) +
  geom_point(size=4)

gss_cat %>% 
  drop_na(tvhours) %>% 
  group_by(relig) %>% 
  summarise(mean_tv = mean(tvhours)) %>% 
  mutate(relig = fct_reorder(relig, mean_tv)) %>% 
  ggplot(aes(mean_tv, relig)) +
  geom_point(size=4)

## By the value of the factor ####
gss_cat %>% 
  ggplot(aes(marital))+
  geom_bar()

gss_cat %>% 
  mutate(marital = fct_infreq(marital)) %>% 
  count(marital)

gss_cat %>% 
  mutate(marital = fct_infreq(marital)) %>% 
  mutate(marital = fct_rev(marital)) %>% 
  count(marital)

gss_cat %>% 
  mutate(marital = marital %>% fct_infreq() %>% fct_rev()) %>% 
  ggplot(aes(marital)) +
  geom_bar()

## By lumping ####
gss_cat %>% 
  count(relig, sort = T)

gss_cat %>% 
  mutate(relig = fct_lump(relig, n = 2)) %>% 
  count(relig)

## By reversing the order ####

gss_cat %>% 
  drop_na(age) %>% 
  filter(rincome != "Not applicable") %>% 
  group_by(rincome) %>% 
  summarise(mean_age = mean(age)) %>% 
  ggplot(aes(mean_age, rincome)) + 
  geom_point(size = 4)

gss_cat %>% 
  drop_na(age) %>% 
  filter(rincome != "Not applicable") %>% 
  group_by(rincome) %>% 
  summarise(mean_age = mean(age)) %>% 
  mutate(rincome = fct_rev(rincome)) %>% 
  ggplot(aes(mean_age, rincome)) + 
  geom_point(size = 4)

## By recoding ####
gss_cat %>% 
  mutate(partyid = fct_recode(partyid,
                              "Republican, strong" = "Strong republican",
                              "Republican, weak" = "Not str republican",
                              )) %>% 
  pull(partyid) %>% 
  levels()


## By collapasing ####
gss_cat %>% 
  mutate(partyid = fct_collapse(partyid,
                                other = c('No answer', "Don't know", "Other party"),
                                rep = c("Strong republican", "Not str republican"),
                                ind = c("Ind,near rep", "Independent", "Ind,near dem"),
                                dem = c("Not str democrat", "Strong democrat")
                                )) %>% 
  count(partyid)



# Separate and unite ####
library(tidyverse)
library(gapminder)

# Separate year into century and year
gapminder1 <- gapminder %>% 
                separate(col = year,
                         into = c("century", "year"),
                         sep = 2) %>% # can also use characters
                view()

gapminder1 %>% 
  unite(col = date, 
        century, year,
        sep = "") %>% 
  view()


