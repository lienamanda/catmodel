#' Values generated to sum to user-defined number
#'
#' @param N integer that the generated values will sum to
#' @param M integer that is number of values that will sum to N
#'
#' @return a double that contains values that sum to a defined number
#'
#' @examples
#' randomSum(5, 2)
#' randomSum(1.2, 5)
#'
#' \dontrun{
#' randomSum(4,2.2)
#' }
randomSum <- function(N, M) {
  if(M %% 1 != 0) {stop("M must be an integer")}
  else {
    vec <- rnorm(N, M/N, 0.1)
    vec / sum(vec) * M
  }
}

# determines whether changing the sensitivity or weight parameter exceeds our defined limits
exceeds <- function(increment, type, x, topLimit) { # type = weight or sensitivity (0)
  if(topLimit < 1) {stop("topLimit cannot be less than 1")
  } else {
    error <- FALSE
    if (increment == 0)  { # decrease
      if ((type == 0 && x < topLimit) | (type > 0 && x < 0.01)) {error <- TRUE}
    } else { # increment = 1 = increase
      if ( (type == 0 && x > 1) | (type > 0 && x > 0.99)) {error <- TRUE}
    }
    return(error)
  }
}

#' Optimized fitting of experimental categorization data to prototype model.
#'
#' This function uses a hill-climbing algorithm to adjust the parameters and
#' minimize the fitting error to report the best-fitting parameters, the fit
#' index, and the prediction of the prototype model.
#'
#' @param protoFile .txt file containing two rows defining the two prototypes.
#'   The number of 0s and 1s (separated by spaces) in each row correspond to the
#'   number of weights/traits in the study.
#' @param obsFile .txt file containing the observed data that will be fitted.
#' Each row corresponds to a different subject, condition, or trial, etc. The
#' number of elements in a row will be equal to the number of columns in the
#' stimulations file.
#' @param stmFile .txt file containing all of the exemplars in the study. These
#'   are typically ordered with dimensions as columns and stimuli as rows. The
#' first "n" rows are items that belong to Category A, then the next "n" rows
#' are itmes that belong to Category B. Additional transfer stimuli (no assigned
#' category) would follow.
#' @param upperBound numeric value (must be greater than 1) that describes
#'   highest value the sensitivity parameter can take (14 is recommended)
#'
#' @return output is creation of 4 .txt files (predicted probabilities, fitted
#'   sensitivities, fitted weight parameters, root mean squared deviations
#'   (RMSD))
#'
#' @examples
#' sum("6dProt.txt", "6dObs.txt", "6dStm.txt", 14)
#' sum("prototypes.txt", "observations.txt", "stimulations.txt", 14.5)
#'
#' \dontrun{
#' sum(file1, file2, file3, 0)
#' }
protoModel <- function(protoFile, obsFile, stmFile, upperBound) {

  prototypes<-read.table(protoFile)
  observations<-read.table(obsFile)
  stimulations<-read.table(stmFile)

  # determines number of subjects/trials, dimensions (traits)
  numTraits <- dim(prototypes)[2]
  numStm <- dim(stimulations)[1]
  numObs <- dim(observations)[1]

  # assumes first and second lines are category A and B prototypes, respectively
  protoA <- array(unlist(prototypes[1, numTraits]), dim = numTraits)
  protoB <- array(unlist(prototypes[2, numTraits]), dim = numTraits)

  # insert experimental data into matrix
  stmData <- matrix(unlist(stimulations[1:numStm, 1:numTraits]), numStm, numTraits)
  obsData <- matrix(unlist(observations[1:numObs, 1:numStm]), numObs, numStm)

  # produces random inital weights
  makeWeights <- TRUE

  while (makeWeights) {
    makeWeights <- FALSE
    weightTable <- matrix(randomSum(numTraits, 1.0), numObs , numTraits, byrow = TRUE)
    for (row in seq(numObs)) {
      if (length(which(weightTable[row, 1:numTraits] < 0)) != 0) {makeWeights <- TRUE}
    }
  }

  # initializes list of RMSDs at value of 1000 and sensitivites to a random value (b/w 0-10)
  RMSD <- array(1000, dim = numObs)
  sens <- array(sample(c(1:2), numObs, replace = TRUE), dim = numObs) # can change starting sens here
  prob <- matrix(NA, numObs, numStm)
  whatChange <- array(NA, dim = numObs) # will later determine what change (weight/sensitivit) will be made to each trial

  # the following determines and applies changes to be made to weight or sensitivity in order to minimize RMSD
  # once counter reaches value (user should choose), the program stops trying to minimize RMSD
  # each weight and sensitivity has an equal chance of being altered

  counter <- 0
  while (counter < 1000) { # <-- user chooses when to stop

    newWeights = weightTable
    newSens = sens
    whatChange<-sample(0:numTraits, numObs, replace = TRUE) # randomly determine what change to make for each trial

    for (currentTrial in seq(numObs)) { # for each trial/subject, apply a change and calculate the RMSD
      direction<-sample(c(0:1), 1) # determines whether to increase(1)/decrease(0)
      if (whatChange[currentTrial] == 0) { # if zero, change sensitivitity
        if (!exceeds(direction, whatChange[currentTrial], sens[currentTrial], upperBound)) {
          if (direction == 0) {
            newSens[currentTrial] <- sens[currentTrial] - 0.01
          } else {
            newSens[currentTrial] <- sens[currentTrial] + 0.01
            }
        }
      } else {

        # CHANGE WEIGHT:
        # if the change in weight surpasses limits, terminate the loop and move onto next trial
        # if weight change is possible, make change and randomly choose another weight to compensate with

        if (!exceeds(direction, whatChange[currentTrial], newWeights[currentTrial, whatChange[currentTrial]], upperBound)) {
          if (direction == 0) {
            newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] - 0.01
          } else {
            newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] + 0.01
            }

          secondCounter <- 0
          while (secondCounter < 20 && secondCounter >= 0) { # looks for and applies change to weight randomly chosen to compensate with
            secondCounter <- secondCounter + 1
            changeTo <- sample(seq(numTraits)[!seq(numTraits) %in% whatChange[currentTrial]], 1) # which other weight will be used to compensate
            swapDirection <- seq(0,1)[!seq(0,1) %in% direction]
            if (!exceeds(swapDirection, changeTo, newWeights[currentTrial,changeTo], upperBound)) { # ensures that the compensating weight change doesn't surpass limits
              if (swapDirection == 1) {
                newWeights[currentTrial, changeTo] <- weightTable[currentTrial, changeTo] + 0.01
                secondCounter <- (-1)
              } else {
                newWeights[currentTrial, changeTo] <- weightTable[currentTrial, changeTo] - 0.01
                secondCounter <- (-1)
                }
            }
          }

          if (secondCounter == 20) {
            if (direction == 0) {
              newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] + 0.01
            } else {
                newWeights[currentTrial,whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] - 0.01
              }
        }
      }
    }

    # finds the weighted differences between both prototypes to experimental data and stores
    testRMSD <- array(dim = numObs) # calculated RMSDs to compare to pre-existing RMSDs
    predictedProb <- matrix(NA, numObs, numStm) # predicted probabilities of choosing category A based on model and new parameters

    # calculates distance, similarity, then probability and stores it in a test array (predictedProb)
    # also calculates RMSD and stores in other test array (testRMSD)
    for (eachTrial in seq(numObs)) {

      differencesA = matrix(protoA, numStm, numTraits, byrow = TRUE) - stmData
      differencesB = matrix( protoB, numStm, numTraits, byrow = TRUE ) - stmData
      weightParams = matrix(newWeights[eachTrial, 1:numTraits], numStm, numTraits, byrow = TRUE)

      distancesA <- abs(differencesA) * weightParams
      distancesB <- abs(differencesB) * weightParams

      # sums the distances of each trait, calculates the similiarity, and inserts the results in the first column
      for (row in seq(numStm)) {

        totalDistanceA = sum(distancesA[row, 1:numTraits])
        totalDistanceB = sum(distancesB[row, 1:numTraits])

        distancesA[row, 1] <- exp(-newSens[eachTrial] * totalDistanceA) # similarity of the stimuli to prototype A
        distancesB[row, 1] <- exp(-newSens[eachTrial] * totalDistanceB) # simliarity of stimuli to B
      }

      # calculate and store probability of classifying stimulus as category A, and then RMSD

      finalDistanceA = distancesA[1:numStm, 1]
      finalDistanceB = distancesB[1:numStm, 1]

      predictedProb[eachTrial, 1:numStm] <- finalDistanceA / (finalDistanceA + finalDistanceB)
      testRMSD[eachTrial] <- sqrt(mean((obsData[eachTrial, 1:numStm] - predictedProb[eachTrial, 1:numStm])^2))

    }

    # if any newly calculated RMSD is less than what was previously stored, replace it and subtract 1 from the counter
    # otherwise, add 1 to the counter

    if (length(which(testRMSD < RMSD)) == 0) { counter <- counter + 1 }
    else {
      for (eachChange in which(testRMSD < RMSD)) {
        weightTable[eachChange, 1:numTraits] <- newWeights[eachChange, 1:numTraits]
        sens[eachChange] <- newSens[eachChange]
        prob[eachChange, 1:numStm] <- predictedProb[eachChange, 1:numStm]
      }

      RMSD[which(testRMSD < RMSD)] <- testRMSD[which(testRMSD < RMSD)]
      counter <- 0
    }
  }

  }

  saveTables(sens, weightTable, RMSD, prob)

}

exemModel <- function(exemFile, obsFile, stmFile, upperBound) {

  exemplars <- read.table(exemFile)
  observations <- read.table(obsFile)
  stimulations <- read.table(stmFile)

  # determines number of subjects/trials, dimensions (traits)
  numTraits <- dim(exemplars)[2]
  numStm <- dim(stimulations)[1]
  numObs <- dim(observations)[1]
  numExemplars <- dim(exemplars)[1]

  # assumes that the first half of training set is Category A and the rest is Category B
  # also assumes that there is an equal amt of A and B trainer sets
  exemA <- matrix(unlist(exemplars[1:(numExemplars / 2 ), 1:numTraits]), numExemplars / 2, numTraits)
  exemB <- matrix(unlist(exemplars[(numExemplars / 2 + 1):numExemplars, 1:numTraits]), numExemplars / 2, numTraits)
  numExemA <- dim(exemA)[1]
  numExemB <- dim(exemB)[1]

  # insert experimental data into matrix
  stmData <- matrix(unlist(stimulations[1:numStm, 1:numTraits]), numStm, numTraits)
  obsData <- matrix(unlist(observations[1:numObs, 1:numStm]), numObs, numStm)

  # produces random inital weights
  makeWeights <- TRUE

  while (makeWeights) {
    makeWeights <- FALSE
    weightTable <- matrix(randomSum(numTraits, 1.0), numObs , numTraits, byrow = TRUE)
    for (row in seq(numObs)) {
      if (length(which(weightTable[row, 1:numTraits] < 0)) != 0) {makeWeights <- TRUE}
    }
  }

  # initializes list of RMSDs at value of 1000 and sensitivites to a random value (b/w 0-10)
  RMSD <- array(1000, dim = numObs)
  sens <- array(sample(c(1:2), numObs, replace = TRUE), dim = numObs)
  prob <- matrix(NA, numObs, numStm)
  whatChange <- array(NA, dim = numObs) # will later determine what change (weight/sensitivit) will be made to each trial

  # the following determines and applies changes to be made to weight or sensitivity in order to minimize RMSD
  # once counter reaches value (user should choose), the program stops trying to minimize RMSD
  # each weight and sensitivity has an equal chance of being altered

  counter <- 0
  while (counter < 1000) { # <-- user chooses when to stop

    newWeights = weightTable
    newSens = sens
    whatChange<-sample(0:numTraits, numObs, replace = TRUE) # randomly determine what change to make for each trial

    for (currentTrial in seq(numObs)) { # for each trial/subject, apply a change and calculate the RMSD
      direction<-sample(c(0:1), 1) # determines whether to increase(1)/decrease(0)
      if (whatChange[currentTrial] == 0) { # if zero, change sensitivitity
        if (!exceeds(direction, whatChange[currentTrial], sens[currentTrial], upperBound)) {
          if (direction == 0) {
            newSens[currentTrial] <- sens[currentTrial] - 0.01
          } else {
            newSens[currentTrial] <- sens[currentTrial] + 0.01
            }
        }
      } else {

        # CHANGE WEIGHT:
        # if the change in weight surpasses limits, terminate the loop and move onto next trial
        # if weight change is possible, make change and randomly choose another weight to compensate with

        if (!exceeds(direction, whatChange[currentTrial], newWeights[currentTrial, whatChange[currentTrial]], upperBound)) {
          if (direction == 0) {
            newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[ currentTrial, whatChange[currentTrial]] - 0.01
          } else {
            newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[ currentTrial, whatChange[currentTrial]] + 0.01
            }

          secondCounter <- 0
          while (secondCounter < 20 && secondCounter >= 0) { # looks for and applies change to weight randomly chosen to compensate with
            secondCounter <- secondCounter + 1
            changeTo <- sample(seq(numTraits)[!seq(numTraits) %in% whatChange[currentTrial]], 1) # which other weight will be used to compensate
            swapDirection <- seq(0,1)[!seq(0,1) %in% direction]
            if (!exceeds(swapDirection, changeTo, newWeights[currentTrial,changeTo], upperBound)) { # ensures that the compensating weight change doesn't surpass limits
              if (swapDirection == 1) {
                newWeights[currentTrial, changeTo] <- weightTable[currentTrial, changeTo] + 0.01
                secondCounter <- (-1)
              } else {
                newWeights[currentTrial, changeTo] <- weightTable[currentTrial, changeTo] - 0.01
                secondCounter <- (-1)
                }
            }
          }

          if (secondCounter == 20) {
            if (direction == 0) {
              newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] + 0.01
            } else {
                newWeights[currentTrial,whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] - 0.01
              }
          }
        }
      }
    }

    # finds the weighted differences between both prototypes to experimental data and stores
    testRMSD <- array(dim = numObs) # calculated RMSDs to compare to pre-existing RMSDs
    predictedProb <- matrix(NA, numObs, numStm) # predicted probabilities of choosing category A based on model and new parameters

    # calculates distance, similarity, then probability and stores it in a test array (predictedProb)
    # also calculates RMSD and stores in other test array (testRMSD)
    for (eachTrial in seq(numObs)) {

      differencesA = matrix(protoA, numStm, numTraits, byrow = TRUE) - stmData
      differencesB = matrix(protoB, numStm, numTraits, byrow = TRUE) - stmData
      weightParams = matrix(newWeights[eachTrial, 1:numTraits], numStm, numTraits, byrow = TRUE)

      distancesA <- abs(differencesA) * weightParams
      distancesB <- abs(differencesB) * weightParams

      # sums the distances of each trait, calculates the similiarity, and inserts the results in the first column
      for ( row in seq(numStm) ) {

        totalDistancesA = sum(distancesA[row, 1:numTraits])
        totalDistancesB = sum(distancesB[row, 1:numTraits])

        distancesA[row, 1] <- exp(-newSens[eachTrial] * totalDistancesA) # similarity of the stimuli to prototype A
        distancesB[row, 1] <- exp(-newSens[eachTrial] * totalDistancesB) # simliarity of stimuli to B
      }

      # calculate and store probability of classifying stimulus as category A, and then RMSD

      finalDistancesA = distancesA[1:numStm, 1]
      finalDistancesB = distancesB[1:numStm, 1]
      predictedProb[eachTrial, 1:numStm] <-  finalDistancesA / (finalDistancesA + finalDistancesB)
      testRMSD[eachTrial] <- sqrt(mean((obsData[eachTrial, 1:numStm] - predictedProb[eachTrial, 1:numStm])^2))

    }

    # if any newly calculated RMSD is less than what was previously stored, replace it and subtract 1 from the counter
    # otherwise, add 1 to the counter

    if (length(which(testRMSD < RMSD)) == 0) {counter <- counter + 1}
    else {
      for (eachChange in which(testRMSD < RMSD)) {
        weightTable[eachChange, 1:numTraits] <- newWeights[eachChange, 1:numTraits]
        sens[eachChange] <- newSens[eachChange]
        prob[eachChange, 1:numStm] <- predictedProb[eachChange, 1:numStm]
      }

      RMSD[which(testRMSD < RMSD)] <- testRMSD[which(testRMSD < RMSD)]
      counter <- 0
    }
  }

  saveTables(sens, weightTable, RMSD, prob)

}

saveTables <- function(sensitivity, weights, RMSDs, probabilities) {
  {title <- readline( prompt = "File name prefix: ")}

  write.table(sensitivity, file = paste0( title, "Sens.txt" ), col.names = F, row.names = F)
  write.table(weights, file = paste0( title, "Weights.txt" ), col.names = F, row.names = F)
  write.table(RMSDs, file = paste0( title, "RMSD.txt" ), col.names = F, row.names = F)
  write.table(probabilities, file = paste0( title, "Probs.txt" ), col.names = F, row.names = F)
}

graphFunc <- function(observed, predicted, whichTrials) {

  listOfTrials <- as.integer(strsplit(whichTrials, " ")[[1]])
  observations <- read.table(observed)
  predictions <- read.table(predicted)
  numStm <- dim(observations)[2]
  stmNumber <- as.double(c(seq(1, numStm)))

  for (eachTrial in listOfTrials) {

    obsProb <- as.double(c(observations[eachTrial, 1:numStm]))
    predProb <- as.double(c(predictions[eachTrial, 1:numStm]))

    dat <- data.frame(
      x_axis = stmNumber,
      Observed = obsProb,
      Predicted = predProb
    )

    reshaped = reshape2::melt(dat, id = c("x_axis"))

    png(paste0("Subject", eachTrial))
    makePlot <- ggplot2::ggplot(reshaped) +
      ggplot2::geom_line(ggplot2::aes(x = x_axis, y = value, colour = variable)) +
      ggplot2::geom_point(ggplot2::aes(x = x_axis, y = value, colour = variable)) +
      ggplot2::xlab("Stimulation") + ggplot2::ylab('Probability of Category "A"') +
      ggplot2::ggtitle(paste0("Subject ", eachTrial)) +
      ggplot2::theme(legend.title = ggplot2::element_blank())
    print(makePlot)
    dev.off()
  }
}
