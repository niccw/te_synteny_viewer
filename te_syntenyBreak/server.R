library(shiny)
library(promises)
library(future)
library(ggplot2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    # all data are in custom environment (data_env)
    id <- NULL # notification id
    warning_id <- NULL # warning notification id
    rv <- reactiveValues(read = 0, basic_plot_param = 0, basic_msynt_plot = NULL, 
                         res_dt = NULL,
                         te_base_cnt1 = NULL, te_base_cnt2 = NULL,
                         ggplot = NULL)
    
    # update file name vector
    speciesA_files <- reactive({
        validate(need(input$speciesA != "","Select species A"))
        print("after validate") # debug
        input_files_by_name(input$speciesA)})
    speciesB_files <- reactive({
        validate(need(input$speciesB != "","Select species B"))
        input_files_by_name(input$speciesB)})
    
    # plot_param_env
    observeEvent({
        plot_param_env$lim <- input$lim
        plot_param_env$nmin <- input$nmin
        plot_param_env$shared <- input$shared
        plot_param_env$showlines <- input$showlines
    },
    {
        print(rv$basic_plot_param)
        rv$basic_plot_param <- rv$basic_plot_param + 1
        print(rv$basic_plot_param)}
    )
    

    #### load files ####

    observeEvent({input$basic_cal},
        {
        validate(need(input$speciesA != "","Select species A"))
        validate(need(input$speciesB != "","Select species B"))
        if(data_env$sp1 != input$speciesA || data_env$sp2 != input$speciesB || is.null(data_env$d1) || is.null(data_env$d2)){
            id <<- showNotification("Reading files...", duration=NULL)
            load_files(filelsA = speciesA_files(), filelsB = speciesB_files())
            removeNotification(id)
            rv$read <- rv$read +1
        }
    })
    
    

    #### basic msynt plot ####
 
    # calculate res_dt
    observeEvent({
        input$basic_cal
    },
                 {
        validate(need(input$speciesA != "","Select species A"))
        validate(need(input$speciesB != "","Select species B"))
        req(data_env$d1,data_env$d2)
        id <<- showNotification("Calculating basic msynt plot", duration=NULL)
        if(data_env$sp1 == input$speciesA && data_env$sp2 == input$speciesB && !is.null(data_env$d1_raw) && !is.null(data_env$d2_raw)){
            rv$basic_msynt_plot <- msynt_matrix_plot(data_env$d1_raw,data_env$d2_raw,input$speciesA,input$speciesB,reactive_env=rv,d_env=data_env)
        }else{
            rv$basic_msynt_plot <- msynt_matrix_plot(data_env$d1,data_env$d2,input$speciesA,input$speciesB,reactive_env=rv,d_env=data_env)
        }
        removeNotification(id)
        id <<- showNotification("Basic msynt dt loaded. Now can plot TE :)", duration=NULL, type="message")
    })
    
    output$basic_plot <- renderPlot({
        validate(need(input$speciesA != "","Select species A"))
        validate(need(input$speciesB != "","Select species B"))
        req(rv$basic_msynt_plot)
        print("inside render plot")
        print(rv$basic_msynt_plot)
    })
    
    output$basic_plot_dl <- downloadHandler(
        filename = function(){
            paste0(input$speciesA,"_",input$speciesB,"_basic_msynt.png")
        },
        content = function(con) {
            png(filename = con,width = input$plot_width,height = input$plot_height, units = "mm", res=96)
            print(rv$basic_msynt_plot)
            dev.off()
        }
    )
    
    # te-msynt
    # TE classes list
    #shinyWidgets::pickerInput("te_classes", label = "TE class",choices = LETTERS ,multiple = TRUE, options = list(`actions-box` = TRUE))
    output$teclass_picker <- renderUI({
        if(input$te_species == "x-axis"){
            print("x-axis")
            te_level <- data_env$te_gff1$type %>% levels() %>% sort()
        }else{
            print("y-axis")
            te_level <- data_env$te_gff2$type %>% levels() %>% sort()
        }
        shinyWidgets::pickerInput("te_classes", label = "TE class",choices = te_level ,multiple = TRUE, options = list(`actions-box` = TRUE))
    })
    
    # count TE and update res_dt
    observeEvent(
        {#obs expression
            validate(
                need(input$update != 0, "Click update to start.")
            )
            input$update
            print("start counting TE")
            },
            # react expression
            {
                req(rv$res_dt) # make sure basic msynt matrix is already here
                
                # new progress object in THIS reactive
                progressA <- shiny::Progress$new()
                progressB <- shiny::Progress$new()
                on.exit({
                    progressA$close()
                    progressB$close()
                    })
                
                # create og to transcript map (only used when considering downstream region as well)
                # d1,d2 scaffold already sorted
                data_env$d1 %<>% separate(V2,c("SP","ID"),sep="\\|")
                data_env$d2 %<>% separate(V2,c("SP","ID"),sep="\\|")
                
                # extract the end position of gene
                # TODO: if use partition, have to extract the end positionof last OG
                data_env$gene_gff1 <- data_env$gene_gff1[data_env$gene_gff1$type == "gene"]
                data_env$gene_gff2 <- data_env$gene_gff2[data_env$gene_gff2$type == "gene"]
                gene_end1 <- GenomicRanges::end(data_env$gene_gff1) # IRange start always smaller than end, don't need to worry the order
                gene_end2 <- GenomicRanges::end(data_env$gene_gff2)
                gene_end1 <- data.table(ID=data_env$gene_gff1$Name, gene_end = gene_end1)
                gene_end2 <- data.table(ID=data_env$gene_gff2$Name, gene_end = gene_end2)
                data_env$d1 %<>% dplyr::full_join(gene_end1)
                data_env$d2 %<>% dplyr::full_join(gene_end2)
                # filter gene not in msynt
                data_env$d1 <- data_env$d1[!is.na(data_env$d1$V1),]
                data_env$d2 <- data_env$d2[!is.na(data_env$d2$V1),]
 
                # A map of Base number to a range of genomic Coord
                #
                if(input$upstream == "5'"){
                    upstream_bool <- TRUE
                }else{
                    upstream_bool <- FALSE
                }
                bc_map1 <- base_coord_map(data_env$d1,scaffold_len1,upstream = upstream_bool,frank_size = input$flank_size)
                bc_map2 <- base_coord_map(data_env$d2,scaffold_len2,upstream = upstream_bool,frank_size = input$flank_size)
                
                # TE count matix, row=base+-window size
                # take some time, for 10000 base, ~3 mins
                # te_base_cnt is the most time consuming part
                progressA$set(message = paste(input$speciesA,"progress:"), value = 0)
                progressB$set(message = paste(input$speciesB,"progress:"), value = 0)
                id <<- showNotification("Counting TE in the flanking region...may take several mins", duration=NULL)
                rv$te_base_cnt1 <- te_base_cnt(bc_map1,data_env$te_gff1,progressA)
                rv$te_base_cnt2 <- te_base_cnt(bc_map2,data_env$te_gff2,progressB)
               
                removeNotification(id)
                id <<- showNotification("Finish TE counting :D", duration=NULL, type="message")
        }
    )
    
    observeEvent({
        validate(
            need(input$update_plot != 0, "Click update to start.")
        )
        input$update_plot

    },{
        #
        req(input$te_classes,input$te_species, rv$te_base_cnt1, rv$te_base_cnt2, rv$res_dt)
        validate(
            need(length(input$te_classes) > 0, "Empty te classes.")
        )
        
        # subset count and tidy up table for plot
        tryCatch({
            print(input$te_classes)
            # FIXME when input$te_classes length=1
            te_df_ls <- te_plot_matrix(tecnt1 = rv$te_base_cnt1, tecnt2 = rv$te_base_cnt2,teClasses = input$te_classes)
        }
            ,
            error = function(e){
                te_df_ls <- NULL
                print(e)
            },
            warning = function(e){
                te_df_ls <- NULL
                print(e)
            }
            )
        
        if (!exists("te_df_ls")){
            rv$ggplot <- NULL
        }else{
            te_df1 <- te_df_ls[1][[1]]
            te_df2 <- te_df_ls[2][[1]]
            tmp_res_dt <- full_join(rv$res_dt,te_df1)
            tmp_res_dt <- full_join(tmp_res_dt,te_df2)
            tmp_res_dt$baseA_n[tmp_res_dt$baseA_n==0] <- NA
            tmp_res_dt$baseB_n[tmp_res_dt$baseB_n==0] <- NA
            tmp_res_dt$baseA_bp[tmp_res_dt$baseA_bp==0] <- NA
            tmp_res_dt$baseB_bp[tmp_res_dt$baseB_bp==0] <- NA
            
            # param
            if(input$te_species == "x-axis"){
                intercept <- "baseA"
                plottype <- paste0("baseA_",input$counttype)
            }else{
                intercept <- "baseB"
                plottype <- paste0("baseB_",input$counttype)
            }
            
            
            # ggplot
            ga <- ggplot(tmp_res_dt,aes(x=baseA,
                                        text = paste("scaffoldA:",scaffoldA,
                                                     "scaffoldB:",scaffoldB,
                                                     "baseNameA:",baseNameA,
                                                     "baseNameB:",baseNameB
                                        )))
            
            if(input$te_species == "x-axis"){
                scaffold_p <- ga + geom_vline(aes_string(xintercept=intercept,color = plottype)) +
                    scale_color_gradient(low = "#ffffff", high = "#ff0000",na.value="white")
            }else{
                scaffold_p <- ga + geom_hline(aes_string(yintercept=intercept,color = plottype)) +
                    scale_color_gradient(low = "#ffffff", high = "#2840fa",na.value="white")
            }
            rv$ggplot <- scaffold_p + 
                geom_point(aes(y=baseB),shape = ".",size=0.3) + 
                theme(panel.background=element_blank(),
                      plot.margin = unit(c(1,1,3,3),"line"))
            if(plot_param_env$showlines){
                rv$ggplot <- rv$ggplot +
                    geom_vline(xintercept = as.numeric(data_env$abx),color="black",linetype="dotted",alpha=0.5, size=0.3) +
                    geom_hline(yintercept = as.numeric(data_env$aby),color="black",linetype="dotted",alpha=0.5, size=0.3)
            }
            rv$ggplot <- rv$ggplot +
                labs(colour = paste("TEcnt_ws:", input$flank_size),
                     x = input$speciesA, y = input$speciesB,
                     title =  paste0("lim=",plot_param_env$lim,";nmin=",plot_param_env$nmin,";TEflank:",input$upstream,"_",input$flank_size,";type:",input$counttype),
                     subtitle = paste(input$te_classes,collapse = ";")) 
        }
    })
    
    
    # render ggplot
    output$static_te_plot <- renderPlot({
        req(rv$ggplot)
        print("inside static_te_plot")
        rv$ggplot
    })
    
    ## download ggplot
    output$ggplot_dl <- downloadHandler(
        filename = function(){
            paste0(input$speciesA,"_",input$speciesB,"_te_msynt.png")
        },
        content = function(con) {
            ggsave(con, plot = rv$ggplot, device = "png", width = input$plot_width, height = input$plot_height, units = "mm")
        }
    )
    
    # render plotly
    output$plotly_o <- renderPlotly({
        input$update_plotly
        req(rv$ggplot)
        print("inside plotly_o")
        isolate(ggplotly(rv$ggplot,tooltip = c("text", "x", "y")))
    })
})
