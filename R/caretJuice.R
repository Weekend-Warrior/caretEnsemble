#' @importFrom magrittr '%>%'
#' @importFrom pipeR '%>>%'
#' @export

# blender
blender <- function(x, y, method = 'glmnet', trControl, ...) {
  suppressWarnings(x %>%
    as.data.frame %>%
    onehot::onehot(max_levels = length(y)) %>>%
    (~ onehot_encoding) %>>%
    predict(x %>%
              as.data.frame,
            sparse = T) %>%
    list(x = .) %>%
    c(list(y = y,
           method = method,
           trControl = trControl),
      list(...)) %>%
    do.call('train', .)) -> blender

  blender$onehot <- onehot_encoding
  return(blender)
}

# caretBlender
# Creates a caretList of univariate encodings
#' @export

caretBlender <- function (x, y, method = 'glmnet', trControl, ...) {
  if (is.null(trControl$indexOut)) stop('caretBlender requires both an index and an indexOut in the trainControl.')

  blendered <- purrr::map(x,
                          ~ blender(.x, y,
                                    method, trControl, ...)) %>%
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

    data <- data[, -conflicts]
  }

  predobs <- makePredObsMatrix(blender)

  #Build a caret model
  model <- train(cbind(as.data.frame(predobs$preds),
                       data[predobs$samps$rowIndex, ]),
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
  stopifnot(any(class(model) == 'caretBlender'))

  if (all(names(model) %in% colnames(data))) {
    conflicts <- which(colnames(data) %in% names(model))

    data <- data[, conflicts]
  } else {
    stop('Not all blended variables are present in the data.frame.')
  }

  model %>%
    purrr::map2_df(.y = data,
                   ~ .y %>%
                     as.data.frame %>%
                     predict(.x$onehot, ., sparse = T) %>%
                     predict(.x, ., ...))
}

# predict.caretJuice
# Takes a caretBlender and a data.frame and returns a vector
#' @export
#' @method predict caretJuice

predict.caretJuice <- function(model, data, ...) {
  stopifnot(any(class(model) == 'caretJuice'))

  blendered <- predict(model$models, data)

  if (any(colnames(data) %in% names(model$models))) {
    conflicts <- which(colnames(data) %in% names(model$models))

    data <- data[, -conflicts]
  }

  stacked <- predict(model$ens_model,
                     cbind(blendered, data),
                     ...)
}
