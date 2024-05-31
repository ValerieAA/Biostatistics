#######################################################
# BIOSCI 220: QUANTITATIVE BIOLOGY
# WEEK 12 LAB
#######################################################

# A hashtag at the beginning of a line means it is a comment and you do not need to run these lines.
# BUT these contain vital information which you should read to help with your understanding of the course and the lab
# The best way to run the code in this lab is to just run line by line. You can either click the run button or control/command + enter. Please avoid highlighting the lines as you may miss a bracket and the code will not work!

#This week we will be making a series of models which will predict the R0 values over different times of the first NZ covid peak. First we have to make the functions and then you will figure out which model is the best for explaining the real data based on likelihood and AIC. 

#######################################################
# LOAD OUR FUNCTIONS FOR SIR SIMULATION AND INFERENCE
#######################################################

library(deSolve)
library(tidyverse)

likelihood_single_point <- function(data.point, model.point, log=TRUE)
	{
	# Get the probability density of the count, assuming the observed count is a poisson process with an average value matching the true process. This allows for variation, e.g. if a case is reported days after it occurs.
	likelihood = suppressWarnings(dpois(x=data.point, lambda=model.point, log=log))
	# Return the likelihood value
	return(likelihood)
	} # END of likelihood_single_point

likelihood_time_series <- function(data.points, model.points, log=TRUE)
	{
	# Error check - return a stop() message if check fails
	if (length(data.points) != length(model.points))
		{ stop("ERROR in likelihood_time_series(): the lengths of data.points and model.points must match.") }
	# Calculate the likelihoods of the data.points
	likelihoods = mapply(FUN=likelihood_single_point, data.point=data.points, model.point=model.points, MoreArgs=list(log=log))
	# If you calculated log-likelihoods, then sum them
	# If you calculated raw likelihoods, then multiply them
	if (log == TRUE)
		{
		total_likelihood = sum(likelihoods)
		} else {
		total_likelihood = prod(likelihoods)  # prod=take the product
		} # END if/else statement: if (log == TRUE)
	# Return that value
	return (total_likelihood)
	} # END of likelihood_time_series

#######################################################
# Functions for drawing an SIR curve from given parameters
#######################################################
# SIR_ode -- this function defines the ordinary differential equations of the SIR model

SIR_ode <- function(time, state, parameters)
	{
	# The above inputs are:
	# time = a time point
	# state = an R list of the 3 state values, i.e. the population sizes of S, I, and R
	# parameters = an R list of the parameters R0 and D_inf
	# Convert the input parameters "R0" and "D_inf" 
	# into the rates "beta" and "nu"
	beta <- parameters[["R0"]] / parameters[["D_inf"]]
	nu <- 1 / parameters[["D_inf"]]
	# Extract the current values of the states S, I, and R
	S <- state[["S"]]
	I <- state[["I"]]
	R <- state[["R"]]
	# Calculate N, the total of S+I+R, 
	# i.e. the population size
	N <- S + I + R
	# Write out the system of 
	# Ordinary Differential Equations (ODEs)
	dS <- -beta * I/N * S 
	dI <- beta * I/N * S - nu * I
	dR <- nu * I
	# Return the current rates of change of S, I, and R:
	return(list(c(dS, dI, dR)))
	} # End of function SIR_ode

# simulate_SIR() -- this runs the simulation, with inputs:
# * theta -- an R list, containing the parameters R0 and D_inf
# * init.state -- the starting counts for S, I, and R individuals
# * times -- the time-points at which to measure the epidemic

simulate_SIR <- function(theta, init.state, times)
	{
	trajectory <- data.frame(ode(y = init.state,
                times = times,
                func = SIR_ode,
                parms = theta,
                method = "ode45"))
	return(trajectory)
	} # End of simulate_SIR()

# simulate_SIR_changes() -- this runs an SIR model in which parameters can change at specific times
# * thetas_states_df, is an R data.frame
# containing the states and parameters at the start, followed by rows describing new parameters at any time points you like.
# * times, are the times to record the epidemic trajectory at

simulate_SIR_changes <- function(thetas_states_df, times)
	{
	# Loop through the rows
	trajectory = NULL
	projected_new_cases = NULL
	projected_recoveries = NULL
	for (rn in 1:nrow(thetas_states_df))
		{
		if (rn == 1)
			{
			init.state = c(S=thetas_states_df$S[rn], I=thetas_states_df$I[rn], R=thetas_states_df$R[rn])
			theta = list(R0=as.numeric(thetas_states_df$R0[rn]), D_inf=as.numeric(thetas_states_df$D_inf[rn]))
			} else {
			last_row = trajectory[nrow(trajectory),]
			init.state = c(S=last_row$S+thetas_states_df$S[rn], I=last_row$I+thetas_states_df$I[rn], R=last_row$R+thetas_states_df$R[rn])
			theta = list(R0=as.numeric(thetas_states_df$R0[rn]), D_inf=as.numeric(thetas_states_df$D_inf[rn]))			
			}
		# Read in the parameters for the initial, versus later, time-bins
		if (rn >= nrow(thetas_states_df))
			{
			start_time = thetas_states_df$time[rn]
			end_time = max(times)
			} else {
			start_time = thetas_states_df$time[rn]
			end_time = thetas_states_df$time[rn+1]
			}
		TF1 = times >= start_time
		TF2 = times <= end_time
		TF = (TF1 + TF2) == 2
		tmp_times = times[TF]
		tmp_trajectory = simulate_SIR(theta=theta, init.state=init.state, times=tmp_times)
		if (rn == 1)
			{
			trajectory = rbind(trajectory, tmp_trajectory)
			} else {
			trajectory_minus_last_state = trajectory[-nrow(trajectory),]
			trajectory = rbind(trajectory_minus_last_state, tmp_trajectory)
			} # if (rn == 1)
		} # END for (rn in 1:nrow(thetas_states_df))
	# Add the projected net new cases
	projected_net_new_cases	= trajectory$I[2:length(trajectory$I)] - trajectory$I[1:(length(trajectory$I)-1)]
	projected_net_new_cases = c(thetas_states_df$I[1], projected_net_new_cases)
	# Subtract the projected recoveries
	projected_recoveries = trajectory$R[2:length(trajectory$R)] - trajectory$R[1:(length(trajectory$R)-1)]
	projected_recoveries = c(0, projected_recoveries)
	# Add back in the projected recoveries, so you get the total number of new cases
	projected_new_cases = projected_net_new_cases + projected_recoveries
	trajectory = cbind(trajectory, projected_net_new_cases, projected_recoveries, projected_new_cases)
	return(trajectory)
	} # End simulate_SIR_changes()

# Function definition for params_to_likelihood_v2
params_to_likelihood_v2 <- function(params, data.points, thetas_states_df, delay_val=0.0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=NULL)
	{
	if (nrow(thetas_states_df) < 2)
		{
		txt = "STOP ERROR in params_to_likelihood_v2(): Your input parameters table, 'thetas_states_df', must have 2 or more rows in it. Check you thetas_states_df input, and try again."
		stop(txt)
		}
	# The times are just the length of the data.points
	times = 1:length(data.points)
	# (1) input the free parameters
	thetas_states_temp = thetas_states_df
	TF = thetas_states_temp == "free"
	thetas_states_temp[TF] = params
	# Ensure everything is numeric
	thetas_states_temp$time = as.numeric(thetas_states_temp$time)
	thetas_states_temp$R0 = as.numeric(thetas_states_temp$R0)
	thetas_states_temp$D_inf = as.numeric(thetas_states_temp$D_inf)
	thetas_states_temp$S = as.numeric(thetas_states_temp$S)
	thetas_states_temp$I = as.numeric(thetas_states_temp$I)
	thetas_states_temp$R = as.numeric(thetas_states_temp$R)
	# You have to round the "time" variable to the nearest day
	thetas_states_temp$time = round(thetas_states_temp$time, digits=0)
	# Error checks, e.g. for "time" values out of order,
	# "time" outside of min/max
	time_inputs_valid = TRUE
	for (i in 2:length(thetas_states_temp$time))
		{
		if (thetas_states_temp$time[i] <= thetas_states_temp$time[i-1])
			{
			time_inputs_valid = FALSE
			break()
			}
		if (thetas_states_temp$time[i-1] < times[1])
			{
			time_inputs_valid = FALSE
			break()			
			}
		if (thetas_states_temp$time[i-1] > times[length(times)])
			{
			time_inputs_valid = FALSE
			break()			
			}
		if (thetas_states_temp$time[i] < times[1])
			{
			time_inputs_valid = FALSE
			break()			
			}
		if (thetas_states_temp$time[i] > times[length(times)])
			{
			time_inputs_valid = FALSE
			break()			
			}
		}
	
	# If any of time inputs are invalid, return an absurdly low lnL
	if (time_inputs_valid== FALSE)
		{
		lnL=-1e100

		# Print, if desired
		if (printvals == TRUE)
			{
			print(thetas_states_temp)
			cat("log-likelihood = ", lnL, "\n", sep="")
			}

		return(lnL)
		}
	
	# Error check
	if ((max(thetas_states_temp$time)+delay_val) > max(times))
		{
		stop("ERROR in params_to_likelihood_v2: there is a 'thetas_states_df$time'+delay_val greater than found in 'times'")
		}
	# Make sure everything in thetas_states_df is numeric, instead of character
	for (i in 1:ncol(thetas_states_temp))
		{
		thetas_states_temp[,i] = as.numeric(thetas_states_temp[,i])
		}
	# Ad2 the delay (allows time for interventions to show up in the detections)
	thetas_states_temp$time[2:length(thetas_states_temp$time)] = thetas_states_temp$time[2:length(thetas_states_temp$time)]
	# Calculate the trajectory given the parameters, input into model.points
	trajectory = simulate_SIR_changes(thetas_states_df=thetas_states_temp, times)
	
	if (use_just_new_cases == TRUE)
		{
		# prob of observed new cases
		model.points = trajectory$projected_new_cases
		lnL1 = likelihood_time_series(data.points, model.points, log=TRUE)
		
		# prob of recoveries
		model.points = trajectory$projected_recoveries
		lnL2 = likelihood_time_series(observed_recoveries, model.points, log=TRUE)
		lnL = lnL1 + lnL2

		} else {
		# The original way -- modeling the number of active cases
		model.points = trajectory$I
		lnL = likelihood_time_series(data.points, model.points, log=TRUE)
		} # END if (use_just_new_cases == TRUE)
	if (is.finite(lnL)== FALSE)
		{
		lnL=-1e100
		}
	# Print, if desired
	if (printvals == TRUE)
		{
		print(thetas_states_temp)
		cat("log-likelihood = ", lnL, "\n", sep="")
		}

	return(lnL)
	} # END params_to_likelihood_v2

#######################################################
# FINISHED LOADING OUR FUNCTIONS FOR SIR SIMULATION AND INFERENCE
#######################################################

#######################################################
# LOAD & PROCESS DATA FROM OURWORLDINDATA
#######################################################

# Function to help plot your country's data
# This function does the data-processing we did in last weeks lab:
# (1) converts "NA" values to 0.0
# (2) adds up the total number of recovered and active cases each day.
calc_active_recovered_cases <- function(d2, add_zeros_back_to=as.Date("2019-12-31"))
	{
	# Subset to your country, replace "NA" values with 0:
	TF = is.na(d2$new_cases)
	d2$new_cases[TF] = 0
	# Keep a count of active cases, and recovered cases
	active_cases = rep(0.0, times=nrow(d2))
	recovered_cases = rep(0.0, times=nrow(d2))
	net_new_cases = rep(0.0, times=nrow(d2))
	num_recovered_cases = rep(0.0, times=nrow(d2))
	for (i in 1:nrow(d2))
		{
		# Keep track of the starting row number
		if (i <= 6)
			{
			startrow = 1
			endrow = i
			recovered_cases[i] = 0
			num_recovered_cases_today = 0
			} else {
			startrow = startrow + 1
			endrow = i
			sum_of_recovered_cases = sum(d2$new_cases[1:(startrow-1)])
			recovered_cases[i] = sum_of_recovered_cases
			num_recovered_cases_today = recovered_cases[i] - recovered_cases[i-1]
			} # END if/else statement
	
		# Add up the last 14 days of active cases
		sum_of_active_cases_14days = sum(d2$new_cases[startrow:endrow])
	
		# Store the result:
		active_cases[i] = sum_of_active_cases_14days
		
		# Calculate net new cases
		net_new_cases[i] = d2$new_cases[endrow] - num_recovered_cases_today
		num_recovered_cases[i] = num_recovered_cases_today
		} # END for (i in 1:nrow(d2))
	
	d3 = cbind(d2, active_cases, recovered_cases, net_new_cases, num_recovered_cases)
	
	# Add 0s back to 2019-12-31 (like the original OurWorldInData had it)
	zero_row = d3[1,]
	zero_row$total_cases = 0
	zero_row$new_cases = 0
	zero_row$active_cases = 0
	zero_row$recovered_cases = 0
	zero_row$net_new_cases = 0
	d3_with_zero_rows_df = d3
	orig_starting_day = as.Date(zero_row$date)
	current_day = orig_starting_day-1
	if (is.na(add_zeros_back_to) == FALSE)
		{
		if (as.Date(d3$date[1]) > as.Date(add_zeros_back_to))
			{
			days_to_add = as.numeric(as.Date(d3$date[1]) - as.Date(add_zeros_back_to))
			for (j in 1:days_to_add)
				{
				zero_row$date = current_day
				d3_with_zero_rows_df = rbind(zero_row, d3_with_zero_rows_df)
				current_day = current_day - 1
				}
			}
		return(d3_with_zero_rows_df)
		}
	return(d3)
	} # END calc_active_recovered_cases <- function(d2)

# Data setup
## download the wk12_data.csv file from CANVAS into what will be your working
## directory for this assignment, then bring the file into RStudio.
## One way you can do this is to use the code below 
## (a box to choose the file from should open and you can select it from there):
wk12_data <- read_csv(file.choose())
## or you can use the instructions provided in the module 1 course guide, chapter 1.4
## (https://biosci220.github.io/BIOSCI220/r-and-rstudio-1.html#dealing-with-data)

country_name <- "New Zealand"
data <- wk12_data %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y"))
TF <- data$location == country_name
d2 <- data[TF, ]
d3 <- calc_active_recovered_cases(d2)

# Plot the data:
ggplot(d3, aes(x = date)) +
  geom_point(aes(y = active_cases), color = "red", shape = "+") +
  geom_point(aes(y = recovered_cases), color = "green3") +
  labs(x = "Date", y = "Active cases")

# To keep it simple, we are cutting off the data at August 15, 2020 (Day #180 of 2020). This is after elimination had been achieved after the First Wave.
# Cutting the data after a date
maximum_date <- as.Date("2020-08-15")
TF <- as.Date(calc_active_recovered_cases(d2)$date) < maximum_date  # TF = TRUE/FALSE result
d3 <- calc_active_recovered_cases(d2)[TF,]  # this takes just the rows where TF=TRUE

# Plot the data:
ggplot(d3, aes(x = date)) +
  geom_point(aes(y = active_cases), color = "red", shape = "+") +
  geom_point(aes(y = recovered_cases), color = "green3") +
  labs(x = "Date", y = "Active cases")

#######################################################
# Seven models of New Zealand's first wave
#######################################################

#######################################################
# Model 1: A 2-regime model, with R0 fixed to 2.7
#######################################################
times <- 1:nrow(d3)      # Number of days since first day
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case

time <- c(1, case1) #So in this model we only have 2 timebins, day 1 to the day of the first case, then the day of the first case forward. Two timebins mean there will be 2 different R0 values, but in all our models we will set the value for days 1 to first case to be zero. 
R0 <- c(0.0, 2.7)    # initial guesses
D_inf <- c(6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0) # putting in the country's population size
I <- c(0, 1) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df==2.7] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)
lnL_result

# Graphing the model model
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- thetas_states_table[TF]
# Convert to numeric
for (i in 1:ncol(thetas_states_ML))
	{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }

# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases) + 0.1 # Adding +0.1 to prevent y-axis error with log(0)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0.1, maxy), xlab="Day", ylab="Number of individuals (log scale)", log="y")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topleft", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM1 (a 2-regime model, R0 fixed to 2.7); max lnL=", round(lnL_result, 2))
title(titletxt)
#What is this plot showing you? Can you figure it out? 


# Save this model's parameters and log-likelihood
thetas_states_ML_model1 <- thetas_states_ML
total_lnL_Model1 <- lnL_result

#######################################################
# Model M2: A 2-regime model (1 free parameter)
#######################################################
times <- 1:nrow(d3)      # Number of days since first day

# Set up regimes:
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case
time <- c(1, case1) #This time we also have 2 timebins, but now we will let the optim function calculate what it thinks the best R0 value would be for day of the first case forward based on the data. This means it is a free parameter. 
R0 <- c(0.0, 3.0)    # initial guesses. these are the values which optim will start from when finding the best R0. 
D_inf <- c(6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0) # putting in the country's population size
I <- c(0, 1) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Make some parameters into "free" parameters, to be inferred. Now we are setting the 2nd R0 value to be free
thetas_states_df$R0[2] <- "free"
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df=="free"] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)


# Running the Maximum Likelihood search with the optim() function.
# LOOK AT THE OUTPUT THAT PRINTS TO SCREEN!!
ML_results <- optim(par=params, fn=params_to_likelihood_v2, data.points=data.points, thetas_states_df=thetas_states_df, delay_val=0, printvals=TRUE, method="L-BFGS-B", lower=0.0, control=list(fnscale=-1), use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Graphing the ML model
# Take the learned parameters from "ML_results", put them
# into a theta_states data.frame for simulation and plotting
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- ML_results$par
for (i in 1:ncol(thetas_states_ML))
	{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }
thetas_states_ML$time[2:length(thetas_states_ML$time)] <- thetas_states_ML$time[2:length(thetas_states_ML$time)]


# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0, maxy), xlab="Day", ylab="Number of individuals")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topright", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM2 (a 2-regime, 1 parameter model), max lnL = ", round(ML_results$value, 2))
title(titletxt)


# Save this model's parameters and log-likelihood
thetas_states_ML_model2 <- thetas_states_ML
total_lnL_Model2 <- ML_results$value

#######################################################
# Model M3: A 3-regime model (3 free parameters)
#######################################################
times <- 1:nrow(d3)      # Number of days since first day
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case

# Set up regimes:
# The first lockdown was on March 26th (Day 87), but it takes several days to see the effect of lockdown. Here we will say it takes 9 days = day 96 will be our "lockdown" day

#Now it is your job to put the appropriate dates in this time bracket.
#For this model we now want 3 times bins, day 1 to case1, case1 to lockdown, and lockdown forward. Enter these 3 times in the brackets below. Replace time1 with the first time etc. It must be in numerical order. 
time <- c(time1, time2, time3) 
R0 <- c(0.0, 3.0, 0.3)    # so you can see now with 1 additional timebin we will have 1 more R0 to calculate. 
D_inf <- c(6.0, 6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0, 0) # putting in the country's population size
I <- c(0, 1, 0) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Make some parameters into "free" parameters, to be inferred
thetas_states_df$I[2] <- "free" #also this time the number of initial infected people is set as a free parameter. Why do you think this is?
thetas_states_df$R0[2] <- "free"
thetas_states_df$R0[3] <- "free"
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df=="free"] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Running the Maximum Likelihood search with the optim() function.
# LOOK AT THE OUTPUT THAT PRINTS TO SCREEN!!
ML_results <- optim(par=params, fn=params_to_likelihood_v2, data.points=data.points, thetas_states_df=thetas_states_df, delay_val=0, printvals=TRUE, method="L-BFGS-B", lower=0.0, control=list(fnscale=-1), use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Graphing the ML model
# Take the learned parameters from "ML_results", put them
# into a theta_states data.frame for simulation and plotting
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- ML_results$par
for (i in 1:ncol(thetas_states_ML))
	{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }
thetas_states_ML$time[2:length(thetas_states_ML$time)] <- thetas_states_ML$time[2:length(thetas_states_ML$time)]

# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0, maxy), xlab="Day", ylab="Number of individuals")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topleft", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM3 (a 3-regime, 3 param model), max lnL = ", round(ML_results$value, 2))
title(titletxt)

# Save this model's parameters and log-likelihood
thetas_states_ML_model3 <- thetas_states_ML
total_lnL_Model3 <- ML_results$value


#######################################################
# Model 4: A 4-regime model (4 free parameters)
#######################################################
times <- 1:nrow(d3)      # Number of days since first day
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case

# Set up regimes:
# The first lockdown was on March 26th (Day 87), but it takes several days to see the effect of lockdown. Here we will say it takes 9 days = 96
# Allowing an uptick in July (in real life, these are imported cases) = 140

#Same here except with another timebin. Read above to figure out what it is! You are adding on an extra day to what you had above. Remember they have to be in numerical order
time <- c(time1, time2, time3, time4) 
R0 <- c(0.0, 3.0, 0.3, 1.5)    # initial guesses
D_inf <- c(6.0, 6.0, 6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0, 0, 0) # putting in the country's population size
I <- c(0, 1, 0, 0) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0, 0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Make some parameters into "free" parameters, to be inferred
thetas_states_df$I[2] <- "free"
thetas_states_df$R0[2] <- "free"
thetas_states_df$R0[3] <- "free"
thetas_states_df$R0[4] <- "free"
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df=="free"] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)


# Running the Maximum Likelihood search with the optim() function.
# LOOK AT THE OUTPUT THAT PRINTS TO SCREEN!!
ML_results <- optim(par=params, fn=params_to_likelihood_v2, data.points=data.points, thetas_states_df=thetas_states_df, delay_val=0, printvals=TRUE, method="L-BFGS-B", lower=0.0, control=list(fnscale=-1), use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Graphing the ML model
# Take the learned parameters from "ML_results", put them
# into a theta_states data.frame for simulation and plotting
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- ML_results$par
for (i in 1:ncol(thetas_states_ML))
	{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }
thetas_states_ML$time[2:length(thetas_states_ML$time)] <- thetas_states_ML$time[2:length(thetas_states_ML$time)]


# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0, maxy), xlab="Day", ylab="Number of individuals")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topleft", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM4 (a 4-regime model, 4 free params) max lnL = ", round(ML_results$value, 2))
title(titletxt)


# Save this model's parameters and log-likelihood
thetas_states_ML_model4 <- thetas_states_ML
total_lnL_Model4 <- ML_results$value

#######################################################
# Model M5: A 5-regime model, 5 free parameters
#######################################################
times <- 1:nrow(d3)      # Number of days since first day
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case

# Set up regimes:
# The first lockdown was on March 26th (Day 87), but it takes several days to see the effect of lockdown. Here we will say it takes 9 days = 96
# Allowing an uptick in July (in real life, these are imported cases) = 140
# Allowing an earlier pre-lockdown slowdown due to social distancing & public health = 88

#one more timebin! (remember to add these in the numerical order to the time object below)

time <- c(time1, time2, time3, time4, time5) 
R0 <- c(0.0, 3.0, 1.5, 0.3, 1.5)    # initial guesses
D_inf <- c(6.0, 6.0, 6.0, 6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0, 0, 0, 0) # putting in the country's population size
I <- c(0, 1, 0, 0, 0) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0, 0, 0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Make some parameters into "free" parameters, to be inferred
thetas_states_df$I[2] <- "free"
thetas_states_df$R0[2] <- "free"
thetas_states_df$R0[3] <- "free"
thetas_states_df$R0[4] <- "free"
thetas_states_df$R0[5] <- "free"
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df=="free"] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)


# Running the Maximum Likelihood search with the optim() function.
# LOOK AT THE OUTPUT THAT PRINTS TO SCREEN!!
ML_results <- optim(par=params, fn=params_to_likelihood_v2, data.points=data.points, thetas_states_df=thetas_states_df, delay_val=0, printvals=TRUE, method="L-BFGS-B", lower=0.0, control=list(fnscale=-1), use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Graphing the ML model
# Take the learned parameters from "ML_results", put them
# into a theta_states data.frame for simulation and plotting
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- ML_results$par
for (i in 1:ncol(thetas_states_ML))
	{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }
thetas_states_ML$time[2:length(thetas_states_ML$time)] <- thetas_states_ML$time[2:length(thetas_states_ML$time)]


# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0, maxy), xlab="Day", ylab="Number of individuals")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topleft", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM5 (a 5-regime, 5 param model) max lnL = ", round(ML_results$value, 2))
title(titletxt)

# Save this model's parameters and log-likelihood
thetas_states_ML_model5 <- thetas_states_ML
total_lnL_Model5 <- ML_results$value


#######################################################
# Model M6: A 6-regime model, 6 free parameters
#######################################################
times <- 1:nrow(d3)      # Number of days since first day
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case

# Set up regimes:
# The first lockdown was on March 26th (Day 87), but it takes several days to see the effect of lockdown. Here we will say it takes 9 days = 96
# Allowing an uptick in July (in real life, these are imported cases) = 140
# Allowing an even earlier pre-lockdown slowdown due to social distancing & public health = 88
# Allowing an earlier pre-lockdown slowdown due to social distancing & public health = 93

# Another one!
time <- c(time1, time2, time3, time4, time5, time6) 
R0 <- c(0.0, 3, 1.5, 1.3, 0.3, 1.5)    # initial guesses
D_inf <- c(6.0, 6.0, 6.0, 6.0, 6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0, 0, 0, 0, 0) # putting in the country's population size
I <- c(0, 1, 0, 0, 0, 0) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0, 0, 0, 0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Make some parameters into "free" parameters, to be inferred
thetas_states_df$I[2] <- "free"
thetas_states_df$R0[2] <- "free"
thetas_states_df$R0[3] <- "free"
thetas_states_df$R0[4] <- "free"
thetas_states_df$R0[5] <- "free"
thetas_states_df$R0[6] <- "free"
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df=="free"] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)


# Running the Maximum Likelihood search with the optim() function.
# LOOK AT THE OUTPUT THAT PRINTS TO SCREEN!!
ML_results <- optim(par=params, fn=params_to_likelihood_v2, data.points=data.points, thetas_states_df=thetas_states_df, delay_val=0, printvals=TRUE, method="L-BFGS-B", lower=0.0, control=list(fnscale=-1), use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Graphing the ML model
# Take the learned parameters from "ML_results", put them
# into a theta_states data.frame for simulation and plotting
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- ML_results$par
for (i in 1:ncol(thetas_states_ML))
	{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }
thetas_states_ML$time[2:length(thetas_states_ML$time)] <- thetas_states_ML$time[2:length(thetas_states_ML$time)]


# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0, maxy), xlab="Day", ylab="Number of individuals")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topleft", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM6 (a 6-regime, 6 param model) max lnL = ", round(ML_results$value, 2))
title(titletxt)

# Save this model's parameters and log-likelihood
thetas_states_ML_model6 <- thetas_states_ML
total_lnL_Model6 <- ML_results$value

#######################################################
# Model M7: A 7-regime model, 7 free parameters
#######################################################
times <- 1:nrow(d3)      # Number of days since first day
case1 <- (1:nrow(d3))[d3$active_cases>0][1] # day of the first detected case

# Set up regimes:
# The first lockdown was on March 26th (Day 87), but it takes several days to see the effect of lockdown. Here we will say it takes 9 days = 96
# Allowing an uptick in July (in real life, these are imported cases) = 140
# Allowing an earlier pre-lockdown slowdown due to social distancing & public health = 88
# Allowing an earlier pre-lockdown slowdown due to social distancing & public health = 93
# Perhaps after the first case and between the lockdown the R0 was even higher as more and more people were getting infected? = 65

# final model. 7 different time sets now

time <- c(time1, time2, time3, time4, time5, time6, time7) 
R0 <- c(0.0, 3, 3.5, 1.5, 1.3, 0.3, 1.5)    # initial guesses
D_inf <- c(6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0) # seems to fit data
S <- c(as.numeric(d3$population[1]), 0, 0, 0, 0, 0, 0) # putting in the country's population size
I <- c(0, 1, 0, 0, 0, 0, 0) # the number of cases on day 1 is a "nuisance parameter"
R <- c(0, 0, 0, 0, 0, 0, 0) # no vaccinations
thetas_states_table <- cbind(time, R0, D_inf, S, I, R)
thetas_states_df <- as.data.frame(thetas_states_table, stringsAsFactors=FALSE)
thetas_states_df

# Make some parameters into "free" parameters, to be inferred
thetas_states_df$I[2] <- "free"
thetas_states_df$R0[2] <- "free"
thetas_states_df$R0[3] <- "free"
thetas_states_df$R0[4] <- "free"
thetas_states_df$R0[5] <- "free"
thetas_states_df$R0[6] <- "free"
thetas_states_df$R0[7] <- "free"
thetas_states_df

# Run the params_to_likelihood_v2() function once, to see your starting likelihood
data_active_cases <- d3$active_cases
data.points <- d3$new_cases
observed_recoveries <- d3$num_recovered_cases
params <- thetas_states_table[thetas_states_df=="free"] # starting parameter values
lnL_result <- params_to_likelihood_v2(params, data.points, thetas_states_df, delay_val=0, printvals=TRUE, use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)


# Running the Maximum Likelihood search with the optim() function.
# LOOK AT THE OUTPUT THAT PRINTS TO SCREEN!!
ML_results <- optim(par=params, fn=params_to_likelihood_v2, data.points=data.points, thetas_states_df=thetas_states_df, delay_val=0, printvals=TRUE, method="L-BFGS-B", lower=0.0, control=list(fnscale=-1), use_just_new_cases=TRUE, observed_recoveries=observed_recoveries)

# Graphing the ML model
# Take the learned parameters from "ML_results", put them
# into a theta_states data.frame for simulation and plotting
thetas_states_ML <- thetas_states_df
TF <- thetas_states_ML == "free"
thetas_states_ML[TF] <- ML_results$par
for (i in 1:ncol(thetas_states_ML))
{ thetas_states_ML[,i] = as.numeric(thetas_states_ML[,i]) }
thetas_states_ML$time[2:length(thetas_states_ML$time)] <- thetas_states_ML$time[2:length(thetas_states_ML$time)]


# Plot the results (note: this plotting is done using baseR, rather than ggplot which we are used to)
trajectory <- simulate_SIR_changes(thetas_states_df=thetas_states_ML, times=times)
maxy <- max(max(trajectory$I), max(data_active_cases))
xvals <- c(trajectory$time, times)
yvals <- c(trajectory$I, data_active_cases)
plot(x=xvals, y=yvals, pch=".", col="white", xlim=c(0, max(trajectory$time)), ylim=c(0, maxy), xlab="Day", ylab="Number of individuals")
lines(x=trajectory$time, y=trajectory$I, lwd=3, col="firebrick2")
points(times, data_active_cases, col="red", pch="+")
legend(x="topleft", legend=c("Active COVID-19 case count", 'ML-fitted projection of "I" (Infected)'), lty=c("blank", "solid"), lwd=c(1,3), pch=c("+", "."), col=c("red","firebrick2"), cex=0.8)

titletxt <- paste0("ML fit, active cases from: ", country_name, "\nM7 (a 7-regime, 7 param model) max lnL = ", round(ML_results$value, 2))
title(titletxt)

# Save this model's parameters and log-likelihood
thetas_states_ML_model7 <- thetas_states_ML
total_lnL_Model7 <- ML_results$value


# Calculate AIC and AIC model weights by completing the Excel Sheet in the assignment instructions:

# Type in the log-likelihoods, and inferred parameters, for your models. For the 
# examples above, this will make R print them to screen.

# Print the log-likelihoods to screen
cat(c(total_lnL_Model1, total_lnL_Model2, total_lnL_Model3, total_lnL_Model4, total_lnL_Model5, total_lnL_Model6, total_lnL_Model7), sep="\n")

# Print the I_ini to screen
cat(c(thetas_states_ML_model1$I[2], thetas_states_ML_model2$I[2], thetas_states_ML_model3$I[2], thetas_states_ML_model4$I[2], thetas_states_ML_model5$I[2], thetas_states_ML_model6$I[2], thetas_states_ML_model7$I[2]), sep="\n")

# Print the R0 to screen
thetas_states_ML_model1$R0
thetas_states_ML_model2$R0
thetas_states_ML_model3$R0
thetas_states_ML_model4$R0
thetas_states_ML_model5$R0
thetas_states_ML_model6$R0
thetas_states_ML_model7$R0


#######################################################
# Congratulations! You have just done basic statistical 
# model comparison.
# To understand the basics of what AIC is about, please
# see the lectures, and read the notes below.
#######################################################




#EXTRA FOR EXPERTS

#######################################################
# Statistical model comparison
#
# One huge advantage of likelihood-based inference is
# that the likelihoods can be used to statistically 
# compare model fits.
#
# One very common method for comparing ML-fitted models 
# is Akaike Information Criterion (AIC). AIC has been 
# used in tens of thousands of research papers to 
# compare multiple model fits. See the lecture slides/
# audio, and this week's reading, for an introduction to 
# AIC.
#
# The very-very short version of AIC is that AIC
# *estimates the relative distance from the truth*.
#
# We never know what the "true model" is in real life,
# since any true biological history that produced observations
# is preposterously complex.  However, we can at least
# compare the *relative* fit of models to data, using AIC.
#
# AIC is a simple function of likelihood:
#
# AIC = -2 * (lnL - nparams)
#
# lnL = the maximized log-likelihood under some model
# nparams = the number of free parameters that you fit with ML
# 
# Basically, the AIC is a form of penalized likelihood: we
# favour models that give higher lnLs, but then we penalize
# models based on the number of free parameters that we fit.
#
# This helps enforce parsimony in our explanations. If we 
# put no penalty on extra free parameters, we could always
# fit models with higher and higher likelihood, until we
# ended up with a free parameter for every data point, at
# which point we fit all the data perfectly, but we have 
# explained nothing at all, because the model just restates
# the dataset.
#
# Some details:
# 
# The "-2" in the AIC formula needs explanation. The "-"
# just converts the log-likelihood (which is always negative,
# if the data are discrete) from negative to positive. This
# makes sense because AIC is a distance, and distances 
# should be positive.
# 
# The "2" part in AIC is actually an arbitrary convention - 
# it would work the same if we used 1, 3, etc. The number "2"
# was chosen because another famous statistical test, the 
# Likelihood Ratio Test, uses 2*(difference in log-likelihoods)
# as the test statistic, often called the "deviance".
#
# The rationale for 1 parameter = 1 unit of log-likelihood
# is complex, but it can be derived from a view of the world 
# where we think that processes are controlled by a series of
# parameters of decreasing importance, such that there are
# a few major controlling variables, and an exponential
# decay in the importance of other controlling variables.
#
# Other rationales can produce other penalties, creating
# other model comparison criteria, such as BIC (Bayesian
# Information Criterion). Don't worry about these for now.
#
# AIC model weights
# 
# AIC is just a number that doesn't mean much by itself.
# AIC values are used to *compare* multiple models. We 
# do this with AIC weights, which assign each model a
# weight, where all of the weights sum to 1.
# 
# For example, if you had these weights:
#
# Model 1 AIC weight = 0.05
# Model 2 AIC weight = 0.9
# Model 3 AIC weight = 0.05
#
# You would say that "Model 2 acquired 90% of the AIC
# weight, and was the best-fitting model in our 
# model set. It fits about 9 times better than the
# other models combined, which is moderately good
# evidence that it is a better fit. Models 1 and 3
# have equal model weight, so they are approximately
# equivalent fits to the data."
#
##############################################
# Here is how to calculate AIC weights:
##############################################
# 
# 1. Calculate AIC values for all the models
#    in your model set.
# 
# 2. Find the AIC-best model. This is the model
#    with the *LOWEST* AIC (i.e., the smallest
#    relative distance from the truth).
#
# 3. For all the models, calculate the AIC difference
#    from the best model. (The AIC difference between
#    the best model and itself is 0.0.).  This
#    number is called "deltaAIC" (in science, "delta" 
#    often means "change" or "difference", so "deltaAIC" 
#    just means "difference in AIC".
# 
# 4. For each model, calculate exp(-0.5 * deltaAIC). These
#    are the relative likelihoods.
#
# 5. Sum the relative likelihoods, and divide each
#    individual relative likelihood by the sum of 
#    relative likelihoods. This is the AIC weight.
#    The AIC weights will sum to 1.0.
#
# 6. If desired, convert these fractions (numbers 
#    between 0 and 1) to percentages, by multiplying
#    by 100, or by clicking the "%" option in Excel.
#
###############################################


###############################################
# Advantages and disadvantages of AIC
###############################################
# The advantages of the AIC framework include:
#
# * Multiple models (more than 2) can be compared easily.
#   Traditional statistical tests, which give you p-values,
#   typically only compare two models at a time, namely
#   a null model and an alternate.
#
# * There is no concept of "p-value cutoffs" or 
#   "statistical signficance" in AIC-based analyses.
#   The proponents of AIC acknowledge that evidence is
#   continuous. If AIC favours a model 1 over model 2
#   60%-40%, this is a tiny bit of evidence that the 
#   model 1 is closer to the truth.  If AIC favours a 
#   model 1 99.999% to 0.001%, then this is strong
#   evidence that model 1 closer to the truth.
#
# * AIC assumes that we never have access to the 
#   complete, true model that produced our data. It
#   acknowledges that the "true" model behind most
#   empirical data would be fantastically complex. So
#   all we are really trying to do is find better 
#   simple approximations of a very complex reality. 
#   AIC is one easy way to measure the *relative* 
#   distance of models from the truth. We never know
#   the "absolute" distance from the truth, because
#   that would require knowing the true model, which
#   we never do.
#
# * Assuming you have maximized lnL values for each
#   model, the formula for AIC, and AIC weights, is very
#   simple.  The calculations can be done by hand, or
#   in Excel or Google Sheets (or R). All you have to 
#   be able to do is count the number of free parameters
#   in each model.
#
# * AIC can be used very broadly for comparing models.
#   This includes not just models where the lnL is 
#   explicitly reported (like those above), but 
#   also in linear regression and other linear models.
#   Here, model fit is more commonly reported in 
#   terms of "deviance" or "RSS" (residual sum
#   of squares, the sum of the squared errors between
#   the line of best fit, and the data).  These are
#   all proportional to -1 * log-likelihood, under
#   the assumption that the errors are normally 
#   distributed. So, the "least squares" method can
#   also be interpreted as a Maximum Likelihood method!
#
#   (Unfortunately, many introductory students only
#    ever learn about least-squares, when Maximum 
#    Likelihood has much broader application to 
#    all sorts of models, not just cases where
#    a normal distribution of errors can be assumed.)
#
# AIC has some limitations as well:
#
# * If you are literally doing a randomized controlled
#   experiment, where you have two hypotheses (models)
#   that exhaust the space of models, namely "no effect"
#   (the null model) and "effect" (the alternative model),
#   then the traditional p-value framework can make a 
#   lot of sense. This is especially the case where
#   further decisions will be based on the result, so 
#   the effect of false-positives and false-negatives
#   needs serious consideration.
#
# * If you *do* have the true model in your set, for
#   example when your data come from a computer simulation
#   where you know the true model, then use a criterion
#   like BIC, which has a stronger penalty against 
#   extra parameters.
#
# * If you have a dataset where the number of data is 
#   small (typically <30 data points), and/or the number
#   of free parameters is large, it is more appropriate
#   to use "AICc", which is sample-size corrected AIC.
#   The formula for this is a little more complex, so
#   we will ignore it here.
#
# * There are all kinds of intense philosophical debates
#   in statistics and philosophy of science about 
#   models, inference, and the best way to measure 
#   support for models/hypotheses. AIC is just one 
#   choice amongst these. I would put it roughly in 
#   the "likelihoodist" school. Other major schools of
#   thought are "frequentism" (which includes p-values
#   and worries about the long-term error rate) and
#   Bayesianism (which explicitly puts prior probability
#   distributions on parameter values, and on models).
#
# * It is easy forget that "all models are wrong", and 
#   get very attached to your best model.  However, this is a 
#   challenge in all of statistics & data analysis.  
#   The famous statement from George Box is
#   "All models are wrong, but some models are useful."
#   Constant critical thinking is required, as I have
#   tried to encourage!
#######################################################



#######################################################
# What is the point of all of this?
#
# In biology, we almost never know the "true model". We
# are just trying to find better approximations of 
# reality.  AIC is one common choice to balance fitting 
# the observations against parsimony.
#
# AIC can be used to address common questions like
# "which model is the best, given the observations I 
#  have", or "out of this set of models, which models
#  are reasonable, given the dataset"?
#
# AIC can be used any time you can get a maximized
# log-likelihood, and a count of the number of free
# parameters, on each model.  It is used so often that
# R has some standard functions to provide them for 
# linear models provided by the lm() function.
# Note: these functions will not work
# outside of standard "base" R functions like lm() and
# glm().
#
# Here is a quick example of that:

data("iris")  # A default R dataset, measurements of iris flowers
plot(iris$Sepal.Length, iris$Petal.Length)

# Run a simple linear model (lm), with Sepal Length
# predicting Petal Length
lm_result <- lm(iris$Petal.Length~Sepal.Length, data=iris)

# The "logLik" function reports the lnL, and the
# number of free parameters (df=3)
# (df means "degrees of freedom", i.e. number of free
#  parameters)
logLik(lm_result)

AIC(lm_result)

# You can see that the AIC formula is being used:
lnL <- as.numeric(logLik(lm_result))
-2*(lnL - 3)
AIC(lm_result)

# 
# MAIN POINT: The skills you have used for using 
# Maximum Likelihood and AIC on SIR models can be
# used anywhere from quite complex process-based 
# models used in epidemiology, ecology, phylogenetics
# etc., all the way to the basic linear models used
# in introductory statistics.
# 
#######################################################





################################################
# Appendix: The connection between Ordinary 
# Least Squares (OLS), i.e. linear regression,
# and Maximum Likelihood.
# 
# This is something I wish I had been taught
# early in my statistics education, so I am 
# mentioning it here.
#
# This Appendix material is "bonus", I will not
# examine you on it for my section of the course.
# 
################################################
# 
# In standard "least squares" analyses
# (OLS, Ordinary Least Squares), where 
# we talk about the line of 
# best fit being the one that minimizes
# the sum of squared residuals (or RSS, residual
# sum of squares), this is just another way
# of talking about maximizing the log-likelihood,
# *if* we make the assumption that the residuals
# (the difference between the predicted line and 
#  the observed response) all follow the same 
# normally distribution (in other words, the 
# residuals are independent, normally distributed
# with constant variance).
# 
# In this situation, the log-likelihood for the
# line of best fit is derived from the log of 
# the probability density function for a normal
# distribution.
# 
# I find that this is rarely explained in R code,
# so I am putting various formulations here.
# 

# Here is the logLik() function code:
getAnywhere("logLik.lm")

# This is the formula from logLik, with weights (w)
# set to 1s (meaning equal weighting of all residuals):
0.5 * - N * (log(2 * pi) + 1 - log(N) + 
        log(sum(lm_result$residuals^2)))


# Here are several simpler formulations that
# give the same result:
# 
# Formula for the log-likelihood of a linear model
# with x predicting y
# Source:
# https://www.stat.cmu.edu/~cshalizi/mreg/15/lectures/06/lecture-06.pdf
# ...page 2, equation (3)
#
# N = number of observations
#
# RSS = residual sum of squares =
# sum of (yi - (b0 + b1*xi))^2
# where
# yi = ith response (observation)
# b0 = intercept
# b1 = slope
# xi = ith predictor
#
# stdev_squared = variance = RSS / N
# stdev = standard deviation = sqrt(stdev_squared) = sqrt(variance)

N <- length(lm_result$residuals)
RSS <- sum(lm_result$residuals^2)
stdev_squared <- RSS/N
stdev <- sqrt(RSS/N)

-N/2 * log(2*pi) - N*log(stdev) - 1/(2*stdev_squared)*RSS

# The formula works the same if you just use RSS
-N/2 * log(2*pi) - N*log(sqrt(RSS/N)) - 1/(2*RSS/N)*RSS

# R's logLik() function has the above formula re-arranged 
# so that RSS appears only once:
0.5 * - N * (log(2 * pi) + 1 - log(N) + 
        log(RSS))

logLik(lm_result)

# These formulations all give the same result!

#######################################################
# References
#
# Cosma Shalizi (2015). "Lecture 6: The Method of 
# Maximum Likelihood for Simple Linear Regression." 
# Lecture 6 from Modern Regression, 36-401, Fall 2015.
# https://www.stat.cmu.edu/~cshalizi/mreg/15/lectures/06/lecture-06.pdf
# 
#######################################################








#######################################################
# References
#
# My SIR code relies heavily on code from the fitR package 
# (Funk et al. 2019).
# 
# Adam, David (2020). "UK has enough intensive care units for 
# coronavirus, expert predicts." New Scientist, March 25, 2020. 
# https://www-newscientist-com.ezproxy.auckland.ac.nz/article/2238578-uk-has-enough-intensive-care-units-for-coronavirus-expert-predicts/
# 
# Anderson, Charles (2020). "New Zealand coronavirus deaths during lockdown could 
# be just 20, modelling suggests." The Guardian. 27 Mar 2020. 
# https://www.theguardian.com/world/2020/mar/27/new-zealand-coronavirus-deaths-during-lockdown-could-be-just-20-modelling-suggests
# 
# Dillane, Tom; Knox, Chris (2020), "Coronavirus Covid 19: New Zealand's 
# best and worst prepared DHBs by elderly population and ICU beds." 
# NZ Herald. April 5, 2020
# https://www.nzherald.co.nz/nz/news/article.cfm?c_id=1&objectid=12321390
# 
# Fefferman, Nina (2020). "The role of applied math in real-time pandemic 
# response: How basic disease models work." Date: 3:30 EDT Tuesday, 
# March 31, 2020.  http://www.nimbios.org/webinars#fefferman 
# 
# The best short introduction to AIC is Franklin et al. (2001), see:
# 
# Franklin, Alan B.; Shenk, Tanya M.; Anderson, David R.; Burnham, Kenneth P. (2001). Statistical model selection: The alternative to hypothesis testing. 
# Chapter 5 of: Modeling in natural resource management: Development, interpretation, and application. Editors: Tanya M. Shenk, Alan B. Franklin. 
# Island Press, Washington, D.C., pp. 75-90. January 2001. 
# 
# PDF at: https://www.researchgate.net/publication/264862335_Statistical_model_selection_The_alternative_to_hypothesis_testing
# Book at: 
# https://books.google.co.nz/books?lr=&id=Uk7rZ7DCvY4C&q=Burhnam#v=snippet&q=Statistical%20model%20selection%3A%20The%20alternative%20to%20hypothesis%20testing&f=false
# https://islandpress.org/books/modeling-natural-resource-management
# 
# Funk, Sebastian (2019). "Introduction to model fitting in R." http://sbfnk.github.io/mfiidd/introduction.html
# 
# Funk, Sebastian; Camacho, Anton; Johnson, Helen; Minter, Amanda; O’Reilly, Kathleen; Davies, Nicholas (2019). "Model fitting and inference for infectious disease dynamics."  Centre for the Mathematical Modelling of Infectious Diseases (CMMID), London School of Hygiene & Tropical Medicine.
# 
# James, Alex; Hendy, Shaun C.; Plank, Michael J.; Steyn, Nicholas (2020). "Suppression and Mitigation 
# Strategies for Control of COVID-19 in New Zealand." 25 March 2020.
# https://cpb-ap-se2.wpmucdn.com/blogs.auckland.ac.nz/dist/d/75/files/2017/01/Supression-and-Mitigation-New-Zealand-TPM-006.pdf
#
# Newton, Kate (2020). "The man modelling NZ's Covid-19 spread 
# from his kitchen table." Radio New Zealand. 27 March 2020.
# https://www.rnz.co.nz/news/in-depth/412744/the-man-modelling-nz-s-covid-19-spread-from-his-kitchen-table
#########################################################