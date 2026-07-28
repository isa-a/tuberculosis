
library(readxl)
library(writexl)
library(dplyr)


df <- read_excel("try.xlsx")


# cat("Dimensions:", dim(df), "\n")
# cat("Column names:", names(df), "\n")
# cat("Sample rows:\n")
# print(head(df, 5))

# ---- 2. DEFINE COUNTRIES AND YEARS ----
# countries <- c(
#   "Afghanistan", "Angola", "Bangladesh", "Benin", "Bhutan", "Botswana",
#   "Burkina Faso", "Burundi", "Cabo Verde", "Cambodia", "Cameroon",
#   "Central African Republic", "Chad", "Comoros", "Congo",
#   "Democratic Republic Of The Congo", "Djibouti", "Equatorial Guinea",
#   "Eritrea", "Eswatini", "Ethiopia", "Gabon", "Gambia", "Ghana",
#   "Greenland", "Guinea", "Guinea-Bissau", "Haiti", "India", "Indonesia",
#   "Ivory Coast", "Kenya", "Kiribati", "Lesotho", "Liberia", "Madagascar",
#   "Malawi", "Mali", "Marshall Islands", "Mauritania", "Mauritius",
#   "Mongolia", "Mozambique", "Myanmar", "Namibia", "Nauru", "Nepal",
#   "Niger", "Nigeria", "North Korea", "Pakistan", "Papua New Guinea",
#   "Peru", "Philippines", "Rwanda", "Sao Tome And Principe", "Senegal",
#   "Seychelles", "Sierra Leone", "Somalia", "South Africa", "South Sudan",
#   "Sudan", "Tanzania", "Thailand", "Timor-Leste", "Togo", "Uganda",
#   "Vietnam", "Zambia", "Zimbabwe"
# )

countries <- c(
  "Albania", "Algeria", "Andorra", "Antigua and Barbuda",
  "Argentina", "Armenia", "Australia", "Austria", "Azerbaijan",
  "Bahamas, The", "Bahrain", "Barbados", "Belarus", "Belgium",
  "Belize", "Bermuda (British)", "Bolivia", "Bosnia and Herzegovina", "Brazil",
  "Brunei", "Bulgaria", "Canada", "Chile", "China",
  "Colombia", "Cook Islands", "Costa Rica", "Croatia", "Cuba",
  "Cyprus", "Czechia", "Ivory Coast", "Denmark", "Dominica",
  "Dominican Republic", "Ecuador", "Egypt", "El Salvador", "Estonia",
  "Fiji", "Finland", "France", "Georgia", "Germany",
  "Greece", "Grenada", "Guatemala", "Guyana",
  "Honduras", "Hungary", "Iceland", "Iran", "Iraq",
  "Ireland", "Israel", "Italy", "Jamaica", "Japan",
  "Jordan", "Kazakhstan", "Kuwait", "Kyrgyzstan", "Latvia",
  "Lebanon", "Libya", "Lithuania", "Luxembourg", "Malaysia",
  "Maldives", "Malta", "Mexico", "Federated States of Micronesia", "Monaco",
  "Montenegro", "Morocco", "Netherlands", "New Zealand", "Nicaragua",
  "Niue", "North Macedonia", "Northern Mariana Islands (United States)", "Norway", "Oman",
  "Palestine", "Panama", "Paraguay", "Poland",
  "Portugal", "Puerto Rico (United States)", "Qatar", "South Korea", "Moldova",
  "Romania", "Russia", "St Kitts and Nevis", "St Lucia", "St Vincent and the Grenadines",
  "Samoa", "San Marino", "Sao Tome and Principe", "Saudi Arabia", "Serbia",
  "Singapore", "Slovakia", "Slovenia", "Solomon Islands", "Spain",
  "Sri Lanka", "Suriname", "Sweden", "Switzerland", "Syria",
  "Taiwan", "Tajikistan", "Tonga", "Trinidad and Tobago",
  "Tunisia", "Turkmenistan", "Tuvalu", "Turkey", "Ukraine",
  "United Arab Emirates", "United States",
  "Uruguay", "Uzbekistan", "Vanuatu", "Venezuela", "Yemen"
)

years <- 2011:2024

# filter info
filtered <- df %>%
  filter(
    Year %in% years,
    Nationality %in% countries,
    `Applicant type` == "All",
    `Case outcome` == "Issued"
  )

# sum dcisions
summary_long <- filtered %>%
  group_by(Nationality, Year) %>%
  summarise(Total_Decisions = sum(Decisions, na.rm = TRUE), .groups = "drop")

# countries as rows, years as columns
summary_wide <- summary_long %>%
  tidyr::pivot_wider(
    names_from  = Year,
    values_from = Total_Decisions,
    values_fill = 0          # fill missing combos with 0
  )


all_countries_df <- data.frame(Nationality = countries)
summary_wide <- all_countries_df %>%
  left_join(summary_wide, by = "Nationality")

# order year cols
year_cols <- as.character(years)
for (yr in year_cols) {
  if (!yr %in% names(summary_wide)) {
    summary_wide[[yr]] <- 0
  }
}
summary_wide[is.na(summary_wide)] <- 0

# reorder cols
summary_wide <- summary_wide %>%
  select(Nationality, all_of(year_cols))

# rename cols
names(summary_wide)[1] <- "Country"


write_xlsx(summary_wide, "decisions_outputLOBURDEN.xlsx")


print(summary_wide)