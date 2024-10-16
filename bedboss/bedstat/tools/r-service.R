processBED = function(path, client, port) {
    message("Processing BED file: ", path)
    
    # Check for shutdown signal
    if (path == "done") {  
        message("Received done signal")
        assign("done", TRUE, envir=.GlobalEnv)
        writeLines("Shutting down server.", con=client) 
        return(0)
    }   
    
    # Check if file exists
    if (!file.exists(path)) {
        message("File not found: ", path)
        writeLines("File not found.", con=client)  
        return(1)
    }

    # Process the file 
    message("File processed successfully")  
    writeLines("File processed successfully", con=client)  
    writeLines("END", con=client)  
    return(1)
}

message("Starting R server")
svSocket::start_socket_server(procfun=processBED)
message ("R server started")
while (!exists("done")) Sys.sleep(1)

message("Shutting down R service")
