# app.R

# Source the global file
source("global.R")

# Source the UI and server files
source("ui.R")
source("server.R")

# Run the app
shinyApp(ui = ui, server = server)
# shiny::runApp()