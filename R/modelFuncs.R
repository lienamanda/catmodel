#' Hill-Climbing Approach to Prototype Model Fitting
#'
#' @param protoFile
#' @param obsFile
#' @param stmFile
#' @param upperBound
#' @param lowerBound
#'
#' @return
#' @export
#'
#' @examples
protoModel <- function(protoFile, obsFile, stmFile, upperBound, lowerBound) {

  prototypes<-read.table(protoFile)
  observations<-read.table(obsFile)
  stimulations<-read.table(stmFile)

  randomSum <- function(N, M, sd = 0.1) {
    vec <- rnorm(N, M/N, sd)
    vec / sum(vec) * M
  }

  # determines whether changing the sensitivity or weight parameter exceeds our defined limits
  exceeds <- function(increment, type, x) { # type = weight or sensitivity (0)
    error <- FALSE
    if (increment == 0)  { # decrease
      if ((type == 0 && x < upperBound) | (type > 0 && x < 0.01)) {error <- TRUE}
    } else { # increment = 1 = increase
      if ( (type == 0 && x > lowerBound) | (type > 0 && x > 0.99)) {error <- TRUE}
    }
    return(error)
  }


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
    weightTable <- matrix( randomSum(numTraits, 1.0), numObs , numTraits, byrow=TRUE)
    for (row in seq(numObs)) {
      if (length(which(weightTable[row, 1:numTraits] < 0)) != 0) { makeWeights<-TRUE }
    }
  }

  # initializes list of RMSDs at value of 1000 and sensitivites to a random value (b/w 0-10)
  RMSD <- array(1000, dim = numObs)
  sens <- array(sample(c(1:2), numObs, replace = TRUE), dim = numObs) # can change starting sens here
  prob <- matrix(NA, numObs, numStm)
  whatChange<-array(NA, dim = numObs) # will later determine what change (weight/sensitivit) will be made to each trial

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
        if (!exceeds(direction, whatChange[currentTrial], sens[currentTrial])) {
          if (direction == 0) { newSens[currentTrial]<-sens[currentTrial] - 0.01}
          else {newSens[currentTrial]<-sens[currentTrial] + 0.01}
        }
      } else {

        # CHANGE WEIGHT:
        # if the change in weight surpasses limits, terminate the loop and move onto next trial
        # if weight change is possible, make change and randomly choose another weight to compensate with

        if (!exceeds(direction, whatChange[currentTrial], newWeights[currentTrial, whatChange[currentTrial]])) {
          if (direction == 0) { newWeights[ currentTrial, whatChange[currentTrial] ] <- newWeights[ currentTrial, whatChange[currentTrial]] - 0.01
          } else {newWeights[ currentTrial, whatChange[currentTrial] ]<-newWeights[ currentTrial, whatChange[currentTrial]] + 0.01}

          secondCounter<-0
          while ( secondCounter < 20 && secondCounter >= 0 ) { # looks for and applies change to weight randomly chosen to compensate with
            secondCounter <- secondCounter + 1
            changeTo <- sample(seq(numTraits)[!seq(numTraits) %in% whatChange[currentTrial]], 1) # which other weight will be used to compensate
            swapDirection <- seq(0,1)[!seq(0,1) %in% direction]
            if (!exceeds(swapDirection, changeTo, newWeights[currentTrial,changeTo])) { # ensures that the compensating weight change doesn't surpass limits
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
            if (direction == 0) {newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] + 0.01}
            else {newWeights[currentTrial,whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] - 0.01}
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
      distancesA <- abs(matrix(protoA, numStm, numTraits, byrow = TRUE) - stmData) * matrix(newWeights[eachTrial, 1:numTraits], numStm, numTraits, byrow = TRUE)
      distancesB<-abs( matrix( protoB, numStm, numTraits, byrow = TRUE ) - stmData ) * matrix(newWeights[ eachTrial, 1:numTraits ], numStm, numTraits, byrow = TRUE)

      # sums the distances of each trait, calculates the similiarity, and inserts the results in the first column
      for ( row in seq(numStm) ) {
        distancesA[row, 1] <- exp(-newSens[eachTrial] * sum(distancesA[row, 1:numTraits])) # similarity of the stimuli to prototype A
        distancesB[row, 1] <- exp(-newSens[eachTrial] * sum(distancesB[row, 1:numTraits])) # simliarity of stimuli to B
      }

      # calculate and store probability of classifying stimulus as category A, and then RMSD
      predictedProb[eachTrial, 1:numStm] <- distancesA[1:numStm, 1] / (distancesA[1:numStm, 1] + distancesB[1:numStm, 1])
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

  printModelFindings(sens, weightTable, RMSD, prob)
  saveTables(sens, weightTable, RMSD, prob)

}

exemModel <- function(exemFile, obsFile, stmFile, upperBound, lowerBound) {

  {
    exemplars <- read.table(exemFile)
    observations <- read.table(obsFile)
    stimulations <- read.table(stmFile)

  }

  randomSum <- function(N, M, sd = 0.1) {
    vec <- rnorm(N, M/N, sd)
    vec / sum(vec) * M
  }

  # determines whether changing the sensitivity or weight parameter exceeds our defined limits
  exceeds <- function(increment, type, x) { # type = weight or sensitivity (0)
    error <- FALSE
    if (increment == 0)  { # decrease
      if ((type == 0 && x < upperBound) | (type > 0 && x < 0.01)) {error <- TRUE}
    } else { # increment = 1 = increase
      if ( (type == 0 && x > lowerBound) | (type > 0 && x > 0.99)) {error <- TRUE}
    }
    return(error)
  }

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
    weightTable <- matrix( randomSum(numTraits, 1.0), numObs , numTraits, byrow=TRUE)
    for (row in seq(numObs)) {
      if (length(which(weightTable[row, 1:numTraits] < 0)) != 0) { makeWeights<-TRUE }
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
        if (!exceeds(direction, whatChange[currentTrial], sens[currentTrial])) {
          if (direction == 0) { newSens[currentTrial]<-sens[currentTrial] - 0.01}
          else {newSens[currentTrial]<-sens[currentTrial] + 0.01}
        }
      } else {

        # CHANGE WEIGHT:
        # if the change in weight surpasses limits, terminate the loop and move onto next trial
        # if weight change is possible, make change and randomly choose another weight to compensate with

        if (!exceeds(direction, whatChange[currentTrial], newWeights[currentTrial, whatChange[currentTrial]])) {
          if (direction == 0) { newWeights[ currentTrial, whatChange[currentTrial] ] <- newWeights[ currentTrial, whatChange[currentTrial]] - 0.01
          } else {newWeights[ currentTrial, whatChange[currentTrial] ]<-newWeights[ currentTrial, whatChange[currentTrial]] + 0.01}

          secondCounter<-0
          while ( secondCounter < 20 && secondCounter >= 0 ) { # looks for and applies change to weight randomly chosen to compensate with
            secondCounter <- secondCounter + 1
            changeTo <- sample(seq(numTraits)[!seq(numTraits) %in% whatChange[currentTrial]], 1) # which other weight will be used to compensate
            swapDirection <- seq(0,1)[!seq(0,1) %in% direction]
            if (!exceeds(swapDirection, changeTo, newWeights[currentTrial,changeTo])) { # ensures that the compensating weight change doesn't surpass limits
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
            if (direction == 0) {newWeights[currentTrial, whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] + 0.01}
            else {newWeights[currentTrial,whatChange[currentTrial]] <- newWeights[currentTrial, whatChange[currentTrial]] - 0.01}
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
      distancesA <- abs(matrix(protoA, numStm, numTraits, byrow = TRUE) - stmData) * matrix(newWeights[eachTrial, 1:numTraits], numStm, numTraits, byrow = TRUE)
      distancesB<-abs( matrix( protoB, numStm, numTraits, byrow = TRUE ) - stmData ) * matrix(newWeights[ eachTrial, 1:numTraits ], numStm, numTraits, byrow = TRUE)

      # sums the distances of each trait, calculates the similiarity, and inserts the results in the first column
      for ( row in seq(numStm) ) {
        distancesA[row, 1] <- exp(-newSens[eachTrial] * sum(distancesA[row, 1:numTraits])) # similarity of the stimuli to prototype A
        distancesB[row, 1] <- exp(-newSens[eachTrial] * sum(distancesB[row, 1:numTraits])) # simliarity of stimuli to B
      }

      # calculate and store probability of classifying stimulus as category A, and then RMSD
      predictedProb[eachTrial, 1:numStm] <- distancesA[1:numStm, 1] / (distancesA[1:numStm, 1] + distancesB[1:numStm, 1])
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

  printModelFindings(sens, weightTable, RMSD, prob)
}

printModelFindings <- function( sensitivity, weights, RMSDs, probabilities ) {
  print("Sensitivities")
  print(sensitivity)
  cat("\n")

  print("Weight Parameters")
  print(weights)
  cat("\n")

  print ("RMSDs")
  print(RMSDs)
  cat("\n")

  print("Probability of Categorizing as A")
  print(probabilities)
  cat("\n")
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
      y_Obs = obsProb,
      y_Pred = predProb
    )

    print(
      ggplot2::ggplot(data = dat, ggplot2::aes(x = x_axis)) +
        ggplot2::geom_line(ggplot2::aes(y = y_Obs), colour = "blue") +
        ggplot2::geom_point(ggplot2::aes(y = y_Obs), colour = "blue") +
        ggplot2::geom_line(ggplot2::aes(y = y_Pred), colour = "red") +
        ggplot2::geom_point(ggplot2::aes(y = y_Pred), colour = "red") +
        ggplot2::xlab("Stimulation") + ggplot2::ylab('Probability of Category "A"') + ggplot2::ggtitle(paste0("Subject ", eachTrial))
    )
  }
}
