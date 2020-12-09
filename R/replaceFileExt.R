replaceFileExt <- function(file.path, ext) {
  sub("\\.[^.]*?$", ext, file.path)
}