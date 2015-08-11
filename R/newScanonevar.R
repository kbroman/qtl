newScanonevar <- function(cross,
                          mean.formula,
                          var.formula,
                          return.effects = FALSE,
                          quiet = TRUE) 
{
  
  validated.input <- ValidateAndFilterInput(cross, mean.formula, var.formula)
  
  scan <- ScanViaDGLM(mean.alt.formula = validated.input$mean.formula,
                      var.alt.formula = validated.input$var.formula,
                      genoprobs = validated.input$genoprobs,
                      mapping.df = validated.input$mapping.df,
                      chr.by.marker = validated.input$chr.by.marker,
                      pos.by.marker = validated.input$pos.by.marker,
                      marker.names = validated.input$marker.names)
  
  return(scan)
}