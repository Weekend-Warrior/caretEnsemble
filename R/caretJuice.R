#' @importFrom magrittr '%>%'
#' @importFrom pipeR '%>>%'
#' @importFrom magrittr '%<>%'

glmnet <- caret::getModelInfo('glmnet')[[2]]

glmnet$predict <- function (modelFit, newdata, submodels = NULL)
{
  if (length(modelFit$obsLevels) < 2) {
    out <- predict(modelFit, newdata, s = modelFit$lambdaOpt)
  }
  else {
    out <- predict(modelFit, newdata, s = modelFit$lambdaOpt,
                   type = "class")
  }
  if (is.matrix(out))
    out <- out[, 1]
  if (!is.null(submodels)) {
    if (length(modelFit$obsLevels) < 2) {
      tmp <- as.list(as.data.frame(predict(modelFit, newdata,
                                           s = submodels$lambda)))
    }
    else {
      tmp <- predict(modelFit, newdata, s = submodels$lambda,
                     type = "class")
      tmp <- if (is.matrix(tmp))
        as.data.frame(tmp, stringsAsFactors = FALSE)
      else as.character(tmp)
      tmp <- as.list(tmp)
    }
    out <- c(list(out), tmp)
  }
  out
}

spectrumString <- caret::getModelInfo('svmSpectrumString')[[1]]

spectrumString$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
  if(any(names(list(...)) == "prob.model") | is.numeric(y))
  {
    out <- ksvm(x = as.list(x), y = y,
                kernel = stringdot,
                kpar = list(type = "spectrum",
                            length = param$length),
                C = param$C, ...)
  } else {
    out <- ksvm(x = as.list(x), y = y,
                kernel = stringdot,
                kpar = list(type = "spectrum",
                            length = param$length),
                C = param$C,
                prob.model = classProbs,
                ...)
  }

  out
}

spectrumString$predict <- function(modelFit, newdata, submodels = NULL) {
  svmPred <- function(obj, x)
  {
    hasPM <- !is.null(unlist(obj@prob.model))
    if(hasPM) {
      pred <- lev(obj)[apply(predict(obj, x, type = "probabilities"),
                             1, which.max)]
    } else pred <- predict(obj, x)
    pred
  }
  out <- try(svmPred(modelFit, as.list(newdata)), silent = TRUE)
  if(is.character(lev(modelFit)))
  {
    if(class(out)[1] == "try-error")
    {
      warning("kernlab class prediction calculations failed; returning NAs")
      out <- rep("", nrow(newdata))
      out[seq(along = out)] <- NA
    }
  } else {
    if(class(out)[1] == "try-error")
    {
      warning("kernlab prediction calculations failed; returning NAs")
      out <- rep(NA, nrow(newdata))
    }
  }
  out
}

# blender
blender <- function(x, y, trControl, ...) {
  UseMethod('blender')
}

#' @method blender factor

blender.factor <- function(x, y, trControl, ...) {
  theDots <- list(...)
  method <- glmnet
  suppressWarnings(x %>%
                     as.data.frame %>%
                     onehot::onehot(max_levels = length(y), addNA = TRUE) %>>%
                     (~ onehot_encoding) %>>%
                     predict(x %>%
                               as.data.frame,
                             sparse = T) %>%
                     list(x = .) %>%
                     c(list(y = y,
                            method = method,
                            trControl = trControl,
                            tuneGrid = expand.grid(alpha = 0,
                                                   lambda = seq(.01, 4.01, .4))),
                       theDots[[1]]) %>%
                     do.call('train', .)) -> blender

  blender$onehot <- onehot_encoding
  return(blender)
}

#' @method blender numeric

blender.numeric <- function(x, y, trControl, ...) {
  theDots <- list(...)
  method <- 'cubist'
  suppressWarnings(x %>%
                     as.data.frame %>%
                     onehot::onehot(max_levels = length(y), addNA = TRUE) %>>%
                     (~ onehot_encoding) %>>%
                     predict(x %>%
                               as.data.frame) %>%
                     list(x = .) %>%
                     c(list(y = y,
                            method = method,
                            trControl = trControl,
                            tuneGrid = expand.grid(committees = 1:3,
                                                   neighbors = 0:2)),
                       theDots[[1]]) %>%
                     do.call('train', .)) -> blender

  blender$onehot <- onehot_encoding
  return(blender)
}

#' @method blender character

blender.character <- function(x, y, trControl, ...) {
  theDots <- list(...)
  method <- spectrumString
  suppressWarnings(x %>%
                     ifelse(is.na(.), "", .) %>%
                     as.data.frame %>%
                     as.matrix %>%
                     list(x = .) %>%
                     c(list(y = y,
                            method = method,
                            trControl = trControl,
                            tuneGrid = expand.grid(length = 2:4,
                                                   C = 10^seq(-3, 3, 2))),
                       theDots[[1]]) %>%
                     do.call('train', .)) -> blender

  return(blender)
}

# caretBlender
# Creates a caretList of univariate encodings
#' @export

caretBlender <- function (x, y, trControl, ...) {
  theDots <- list(...)

  if (is.null(trControl$indexOut)) stop('caretBlender requires both an index and an indexOut in the trainControl.')

  blendered <- purrr::map(x,
                          ~ blender(.x, y, trControl = trControl, theDots)) %>%
    setNames(nm = colnames(x))

  class(blendered) <- c('caretBlender', 'caretList')
  return(blendered)
}

# caretJuice
# Takes a caretBlender and a data.frame, blenders then stacks
#' @export

caretJuice <- function (blender, data, ...) {
  stopifnot(any(class(blender) == 'caretBlender'))

  if (any(colnames(data) %in% names(blender))) {
    conflicts <- which(colnames(data) %in% names(blender))

    if (!all(colnames(data) %in% names(blender))) data <- data[, -conflicts]
  }

  predobs <- makePredObsMatrix(blender)

  if (!all(colnames(data) %in% names(blender))) {
    newx <- cbind(as.data.frame(predobs$preds),
                     data[predobs$samps$rowIndex, ])
  } else {
    newx <- predobs$preds
  }

  #Build a caret model
  model <- train(newx,
                 predobs$obs, ...)

  #Return final model
  out <- list(models = blender,
              ens_model = model,
              error = model$results)
  class(out) <- c('caretJuice', "caretStack")
  return(out)
}

# predict.caretBlender
# Takes a caretBlender and a data.frame and returns a blendered data.frame
#' @export
#' @method predict caretBlender

predict.caretBlender <- function (model, data, ...) {
  theDots <- list(...)

  if (all(names(model) %in% colnames(data))) {
    name_string <- names(model) %>%
      paste(collapse = ', ') %>%
      paste0('c(', ., ')')

    data %<>%
      dplyr::select_(name_string) %>%
      as.data.frame
  } else {
    stop('Not all blended variables are present in the data.frame.')
  }

  model %>%
    purrr::map2_df(.y = data,
                   ~ .y %>%
                     as.data.frame %>%
                     predict(.x$onehot, ., sparse = !is.null(.x$onehot$.$levels)) %>%
                     list(.x, .) %>%
                     c(theDots) %>%
                     do.call('predict', .))
}

# predict.caretJuice
# Takes a caretBlender and a data.frame and returns a vector
#' @export
#' @method predict caretJuice

predict.caretJuice <- function(model, data, ...) {

  blendered <- predict(model$models, data)

  if (any(colnames(data) %in% names(model$models))) {
    conflicts <- which(colnames(data) %in% names(model$models))

    if (!all(colnames(data) %in% names(model$models))) data <- data[, -conflicts]
  }

  if (!all(colnames(data) %in% names(model$models))) {
    newdata <- cbind(blendered, data)
  } else {
    newdata <- blendered
  }

  stacked <- predict(model$ens_model,
                     newdata,
                     ...)
}
