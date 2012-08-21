library(Rook)
s <- Rhttpd$new()
dir.create(file.path(tempdir(),"XML"),showWarnings=FALSE)
s$add( name="ibiscore",
       app=Builder$new(
         Static$new(
           urls = c("/css","/images","/js"),
           root = getwd()
         ),
         Static$new(urls="/XML",root=tempdir()),
         Brewery$new(
           url="/brew",
           root= getwd()
         ),
         Redirect$new("/brew/ibiscore_gmap.rhtml")
       )
)
s$start(quiet=TRUE)
s$browse("ibiscore")