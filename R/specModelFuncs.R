specProtoModel <- function(protoFile, obsFile, stmFile, paramFile) {

  counter <- 1

  prototypes <- read.table(protoFile)
  observations <- read.table(obsFile)
  stimulations <- read.table(stmFile)
  parameters <- read.table(paramFile)   # NOTE** Sensitivities in first column

  # determines number of subjects/trials, dimensions (traits)
  numTraits <- dim(prototypes)[2]
  numStm <- dim(stimulations)[1]
  numObs <- dim(observations)[1]
  numTests <- dim(parameters)[1]

  # assumes that the first half of training set is Category A and the rest is Category B
  # also assumes that there is an equal amt of A and B trainer sets
  protoA <- array(unlist(prototypes[1, numTraits]), dim = numTraits)
  protoB < -array(unlist(prototypes[2, numTraits]), dim = numTraits)

  # insert experimental data into matrix
  stmData <- matrix(unlist(stimulations[1:numStm, 1:numTraits]), numStm, numTraits)
  obsData <- matrix(unlist(observations[1:numObs, 1:numStm]), numObs, numStm)
  parameterData <- matrix(unlist(parameters[1:numTests, 1:(numTraits+1)]), numTests, numTraits+1) # +1 for sensitivity

  totalProb <- matrix(NA, numObs, numStm + 1)

  for (eachTest in seq(numTests)) {
    weights <- array(parameterData[eachTest, 2:(numTraits + 1)], dim = numTraits)
    sens <- parameterData[eachTest, 1]

    predictedProb <- array(dim = numStm) # will contain predicted probabilities of choosing category A based on model and new parameters

    # finds the weighted differences between all prototypes to experimental data and stores
    # by storage and calculations in multi-dimensional arrays, where each 3rd dimension represents a trainer

    distancesA <- abs(matrix(protoA, numStm, numTraits, byrow = TRUE) - stmData) * matrix( weights, numStm, numTraits, byrow = TRUE)
    distancesB <- abs(matrix(protoB, numStm, numTraits, byrow = TRUE) - stmData) * matrix(weights, numStm, numTraits, byrow = TRUE)

    # sums the distances of each trait, calculates the similiarity, and inserts the results in the first column
    for (row in seq(numStm)) {
      distancesA[row, 1] <- exp(-sens * sum(distancesA[row, 1:numTraits])) # similarity of the stimuli to prototype A
      distancesB[row, 1] <- exp(-sens * sum(distancesB[row, 1:numTraits])) # simliarity of stimuli to B    }
    }

    # calculate and store probability of classifying stimulus as category A, and then RMSD
    predictedProb[1:numStm] <- distancesA[1:numStm, 1] / (distancesA[1:numStm, 1] + distancesB[ 1:numStm, 1])
    RMSD <- sqrt(mean((obsData[eachTest, 1:numStm] - predictedProb[1:numStm]) ^ 2))

    totalProb[counter, 1] <- RMSD # RMSD in first column
    totalProb[counter, 2:(numStm + 1)] <- predictedProb[1:numStm]

    counter <- counter + 1
  }

  printSpecModelFindings(totalProb)
}

specExemModel <- function() {

  counter <- 1
  totalProb <- matrix(NA, numObs, numStm + 1)

  exemFile <- readline(prompt = "Enter exemplar text file: ") # text file containing exemplars
  obsFile <- readline(prompt = "Enter observations text file: ")  # text file containing observational data to be fit to
  stmFile <- readline(prompt = "Enter stimulations text file: ")  # text file containing stimulations
  paramFile <- readline(prompt = "Enter parameters text file: ") # text file containing weights and sensitivities

  exemplars <- read.table(exemFile)
  observations <- read.table(obsFile)
  stimulations <- read.table(stmFile)
  parameters <- read.table(paramFile)   # NOTE** Sensitivities in first column

  # determines number of subjects/trials, dimensions (traits)
  numTraits <- dim(exemplars)[2]
  numStm <- dim(stimulations)[1]
  numObs <- dim(observations)[1]
  numExemplars <- dim(exemplars)[1]
  numTests <- dim(parameters)[1]

  # assumes that the first half of training set is Category A and the rest is Category B
  # also assumes that there is an equal amt of A and B trainer sets
  exemA <- matrix(unlist(exemplars[1:(numExemplars / 2), 1:numTraits]), numExemplars / 2, numTraits)
  exemB <- matrix(unlist(exemplars[(numExemplars / 2 + 1):numExemplars, 1:numTraits]), numExemplars / 2, numTraits)
  numExemA <- dim(exemA)[1]
  numExemB <- dim(exemB)[1]

  # insert experimental data into matrix
  stmData <- matrix(unlist(stimulations[1:numStm, 1:numTraits]), numStm, numTraits)
  obsData <- matrix(unlist(observations[ 1:numObs, 1:numStm]), numObs, numStm)
  parameterData <- matrix(unlist(parameters[1:numTests, 1:(numTraits + 1)]), numTests, numTraits + 1) # +1 for sensitivity
  totalProbAndRMSD <- matrix(NA, numObs, numStm + 1)

  for (eachTest in seq(numTests)) {
    weights <- array(parameterData[eachTest, 2:(numTraits + 1)], dim = numTraits)
    sens <- parameterData[eachTest, 1]

    predictedProb <- array(dim = numStm) # will contain predicted probabilities of choosing category A based on model and new parameters

    # finds the weighted differences between all exemplars to experimental data and stores
    # by storage and calculations in multi-dimensional arrays, where each 3rd dimension represents a trainer
    distancesA <- array(dim = c(numStm, numTraits, numExemA))
    distancesB <- array(dim = c(numStm, numTraits, numExemB))

    for (eachTrainer in seq(numExemA) ) { # for each trainer
      distancesA[1:numStm, 1:numTraits, eachTrainer] <- abs(matrix(exemA[eachTrainer,1:numTraits], numStm, numTraits, byrow = TRUE) - stmData) * matrix(weights, numStm, numTraits, byrow = TRUE)
      distancesB[1:numStm, 1:numTraits, eachTrainer] <- abs(matrix(exemB[eachTrainer,1:numTraits], numStm, numTraits, byrow = TRUE) - stmData) * matrix(weights, numStm, numTraits, byrow = TRUE)

      for (row in seq(numStm)) { # calculates similarity of each stimuli to category A and B across all exemplars
        distancesA[row, 1, eachTrainer] <- exp(-sens * sum(distancesA[row, 1:numTraits, eachTrainer]))
        distancesB[row, 1, eachTrainer] <- exp(-sens * sum(distancesB[row, 1:numTraits, eachTrainer]))
      }
    }

    for (eachRow in seq(numStm)) { # calculates means of sums of all exemplars and stores in first column of first matrix
      distancesA[eachRow, 1, 1] <- sum(distancesA[ eachRow, 1, 1:numExemA])
      distancesB[eachRow, 1, 1] <- sum(distancesB[ eachRow, 1, 1:numExemB])
    }

    # calculate and store probability of classifying stimulus as category A, and then RMSD
    predictedProb[1:numStm] <- distancesA[1:numStm, 1, 1] / (distancesA[1:numStm, 1, 1] + distancesB[1:numStm, 1, 1])
    RMSD <- sqrt(mean((obsData[eachTest, 1:numStm] - predictedProb[1:numStm]) ^ 2))

    totalProb[counter, 1] <- RMSD # RMSD in first column
    totalProb[counter, 2:(numStm + 1)] <- predictedProb[1:numStm]

    counter <- counter + 1
  }
  printSpecModelFindings(totalProb)

}

printSpecModelFindings <- function(probAndRMSDTable) {
  print("RMSDs")
  print(probAndRMSDTable[1:dim(probAndRMSDTable)[1], 1])
  cat("\n")

  print("Probabilities of Categorizing as A")
  print(probAndRMSDTable[1:dim(probAndRMSDTable)[1], 2:dim(probAndRMSDTable)[2]])
}
