#' @importFrom magrittr '%>%'
#' @importFrom pipeR '%>>%'
#' @importFrom magrittr '%<>%'

# blender
blender <- function(x, y, trControl, ...) {
  UseMethod('blender')
}

#' @method blender factor

blender.factor <- function(x, y, trControl, ...) {
  theDots <- list(...)
  method <- 'glmnet'
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

  if (all(names(model) %in% colnames(data))) {
    data %<>%
      dplyr::select_(paste0('c(',
                            paste(names(model),
                                  sep = ', '),
                            ')')) %>%
      as.data.frame
  } else {
    stop('Not all blended variables are present in the data.frame.')
  }

  model %>%
    purrr::map2_df(.y = data,
                   ~ .y %>%
                     as.data.frame %>%
                     predict(.x$onehot, ., sparse = !is.null(.x$onehot$.$levels)) %>%
                     predict(.x, ., ...))
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
