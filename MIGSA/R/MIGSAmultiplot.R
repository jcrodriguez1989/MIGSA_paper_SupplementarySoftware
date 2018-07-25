setGeneric(name="multiplot", def=function(plotlist, layout) {
    standardGeneric("multiplot")
})


# Multiple plot function
#
# ggplot objects must be passed in plotlist (as a list of ggplot objects)
# - layout: A matrix specifying the layout.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
#'@importFrom grid grid.layout grid.newpage pushViewport viewport
setMethod(
    f="multiplot",
    signature=c("list", "matrix"),
    definition=function(plotlist, layout) {
        plots <- c(plotlist)
        numPlots = length(plots)
        
        if (numPlots==1) {
            print(plots[[1]])
        } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(nrow(layout), 
                ncol(layout))))
            
            # Make each plot, in the correct location
            for (i in 1:numPlots) {
                # Get the i,j matrix positions of the regions that contain this 
                # subplot
                matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                
                print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                layout.pos.col = matchidx$col))
            }
        }
    }
)
