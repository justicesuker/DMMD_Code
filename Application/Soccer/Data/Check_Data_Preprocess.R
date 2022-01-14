data = read.csv("Application/Soccer/Data/Premier_League_matches.csv")
dim(data)
head(data)
data[1:5,1:5]
data[1,1]
colnames(data)
data_filtered = data[,c("HomeTeam","AwayTeam",
                        "Full.Time.Home.Goals","Full.Time.Away.Goals",
                        "Full.Time.Results",
                        "Half..Time.Home.Goals","Half.Time.Away.Goals",
                        "Home.Shots","Away.Shots",
                        "Home.Shots.on.Target","Away.Shots.on.Target",
                        "Home.Hit.Woodwork","Away.Hit.Woodwork",
                        "Home.Corners","Away.Corners",
                        "Home.Fouls.Committed", "Aways.Fouls.Committed",
                        "Home.Offsides", "Away.Offsides", 
                        "Home.Yellow.Cards", "Away.Yellow.Cards",
                        "Home.Red.Cards", "Away.Yellow.Cards.1")]

# Change the last column name to be "red cards". We believe it's a typo.
colnames(data_filtered) = c("HomeTeam","AwayTeam",
                            "Full.Time.Home.Goals","Full.Time.Away.Goals",
                            "Full.Time.Results",
                            "Half..Time.Home.Goals","Half.Time.Away.Goals",
                            "Home.Shots","Away.Shots",
                            "Home.Shots.on.Target","Away.Shots.on.Target",
                            "Home.Hit.Woodwork","Away.Hit.Woodwork",
                            "Home.Corners","Away.Corners",
                            "Home.Fouls.Committed", "Aways.Fouls.Committed",
                            "Home.Offsides", "Away.Offsides", 
                            "Home.Yellow.Cards", "Away.Yellow.Cards",
                            "Home.Red.Cards", "Away.Red.Cards")
head(data_filtered)
n = dim(data_filtered)[1]
wl_data_filtered = data_filtered[data$"Full.Time.Results"!="D",]
dim(wl_data_filtered)
dim(data)

# The following process will result in the same data that Himanshu gets.
rows_selected = wl_data_filtered$Home.Yellow.Cards == round(wl_data_filtered$Home.Yellow.Cards)&
                wl_data_filtered$Away.Yellow.Cards == round(wl_data_filtered$Away.Yellow.Cards)&
                wl_data_filtered$Home.Red.Cards == round(wl_data_filtered$Home.Red.Cards)&
                wl_data_filtered$Away.Red.Cards == round(wl_data_filtered$Away.Red.Cards)

wl_no_decimal_cards_data = wl_data_filtered[rows_selected,]
dim(wl_no_decimal_cards_data)

home_column = c("Full.Time.Home.Goals","Half..Time.Home.Goals","Home.Shots","Home.Shots.on.Target",
                "Home.Hit.Woodwork","Home.Corners","Home.Fouls.Committed", "Home.Offsides",
                "Home.Yellow.Cards","Home.Red.Cards")

away_column = c("Full.Time.Away.Goals","Half.Time.Away.Goals","Away.Shots","Away.Shots.on.Target",
                "Away.Hit.Woodwork","Away.Corners","Aways.Fouls.Committed","Away.Offsides", 
                "Away.Yellow.Cards","Away.Red.Cards")

home_win_data = wl_no_decimal_cards_data[wl_no_decimal_cards_data$Full.Time.Results == 'H', home_column]
away_win_data = wl_no_decimal_cards_data[wl_no_decimal_cards_data$Full.Time.Results == 'A', away_column]
home_lose_data = wl_no_decimal_cards_data[wl_no_decimal_cards_data$Full.Time.Results == 'A', home_column]
away_lose_data = wl_no_decimal_cards_data[wl_no_decimal_cards_data$Full.Time.Results == 'H', away_column]
# one Home-win match is corresponding to one away-lose match and vice versa. So the following two
# data sets are double-matched.
matched_colnames = c("Full.Time.Goals","Half.Time.Goals","Shots","Shots.on.Target","Hit.Woodwork",
                     "Corners","Fouls.Committed","Offsides", "Yellow.Cards","Red.Cards")
colnames(home_win_data) = matched_colnames
colnames(away_win_data) = matched_colnames
colnames(home_lose_data) = matched_colnames
colnames(away_lose_data) = matched_colnames

win_data = rbind(home_win_data,away_win_data)
lose_data = rbind(away_lose_data,home_lose_data)

write.csv(win_data,"Win_team.csv",row.names = F)
write.csv(lose_data,"Lose_team.csv",row.names = F)