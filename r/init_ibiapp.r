library(Rook)
s <- Rhttpd$new()
dir.create(file.path(tempdir(),"plots"),showWarnings=FALSE)
s$add( name="ibiscore",
       app=Builder$new(
         Static$new(
           urls = c("/css","/images","/js"),
           root = getwd()
         ),
         Static$new(urls="/plots",root=tempdir()),
         Brewery$new(
           url="/brew",
           root= getwd(),
           imagepath=file.path(tempdir(),"plots"),
           imageurl="../plots/"
         ),
         Redirect$new("/brew/ibiscore_gmap.rhtml")
       )
)
s$start(quiet=TRUE)
s$browse("ibiscore")