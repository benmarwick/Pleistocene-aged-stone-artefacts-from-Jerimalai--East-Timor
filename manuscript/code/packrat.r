# Load dependency management system to improve portability and
# reproducibility. This stores all the packages and dependancies 
# in the project folder, rather than my personal library.
# More details: http://rstudio.github.io/packrat/

# just do this once when starting the project
# packrat::init()

# do this periodically to see if I've added any packages that haven't
# been locally stored. Want to do this each time I work on the project.
# packrat::status()

# if the status says I have some packages that aren't yet in the
# packrat. Can run every time I work on the project & add or remove a 
# package
# packrat::snapshot() 

# if the status says it found packages I'm not using, use this to clean them # out
# packrat::clean()

# if this project has been copyied onto a new machine, or if it's been
# updated via version control, then I (or you) need to run this:
# packrat::restore()