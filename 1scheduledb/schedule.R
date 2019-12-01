#use read.csv because the format of the 1st row can maintain as "character"
#which is important that enable dbGetQuery perform smoothly 
raw_data<-read.csv("download.csv",header = T)
library(RSQLite)

# Create a connection to our new database, scheduleDB.db
# you can check that the .db file has been created on your working directory
conn <- dbConnect(RSQLite::SQLite(), "scheduleDB.db")

#Write the raw_data dataset into a table names schedule_data
dbWriteTable(conn, "schedule_data", raw_data,overwrite=T)

#List all the tables available in the database
dbListTables(conn)

# Gather the first 10 rows in the schedule_data table
dbGetQuery(conn, "SELECT * FROM schedule_data LIMIT 10")

#Gather date,week,lesson2 columns' first 10 rows in the schedule_data table
dbGetQuery(conn, "SELECT date,week,lesson2 FROM schedule_data LIMIT 10")

#Get all columns with date 9.3
dbGetQuery(conn,"SELECT * FROM schedule_data WHERE date = 903")

#Get the date,week,lesson2 columns with date 9.3
dbGetQuery(conn,"SELECT date,week,lesson2 FROM schedule_data WHERE date = 903")

#Get the date,lesson2 columns with date pattern including '9'&week='三','四','五'
dbGetQuery(conn,"SELECT date,week,lesson2 FROM schedule_data WHERE date LIKE '9%' AND week IN ('三','四','五')")

# Close the database connection to scheduleDB
dbDisconnect(conn)