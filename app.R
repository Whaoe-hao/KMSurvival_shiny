require(shiny)
require(survival)
require(rsconnect)
library(shinythemes)
ui <-  tagList(
  navbarPage(
    theme = shinythemes::shinytheme("flatly"),
    "Platform for analysis of lipid genes in cancer",
    tabPanel("HOME",
             fluidRow(
               column(12, offset = 0.5,
                      mainPanel(
                        br(),
                        br(),
                        h3("Methods:"), h4("R package Survminer and Survival were used for Kaplan-Meier Survival analysis.  TCGA clinical data were read, and gene expression data and corresponding clinical data were combined according to sample TCGA barcode.  The total survival time of each patient was calculated based on patient survival status (Vital_status),
                      time from diagnosis to death (DAYs_to_death), and time from diagnosis to recent follow-up (days_to_last_FOLLOW_up) in clinical data.  Combined clinical data and gene mutation data.  Based on mutations in specific lipid metabolism genes, the samples were divided into two groups: those that had the mutation and those that did not.  Chi-square
                      test was used to analyze the significance of difference between the two survival curves (corresponding to mutant samples and non-mutant samples), and the results with significant difference (p-value<0.05) were screened out "),
                        br(),
                        br(),
                        h3("About the Data:"), h4("Tumor transcriptome data, sample description data and clinical data were downloaded from TCGA database.  We selected 14 cancer species with both normal samples and tumor sample sizes greater than 19 as our research objects, and their sample information is shown in Figure", strong("(1)"),
                                                  "They were sorted and classified by R language, and the expression data of protein-coding genes were retained to obtain the gene expression data of normal samples and tumor samples.  Combined with the gene expression matrix and sample group information file, normalized gene expression data
                                                were obtained by using the correlation function in R package DESeq2, and differential expression analysis was conducted to screen out differential expression genes.  Using the above methods to analyze and process the sample information, the data results obtained are the data basis of this survival analysis "),
                        
                        br(),
                        br(),
                        img(src="TCGA_sample.png", width=1059, height=644),
                        
                      )
               )
             )
    ),
    tabPanel("Survival analysis",
             fluidRow(
               column(12, offset = 0.5, 
                      sidebarLayout(
                        sidebarPanel(
                          selectInput(inputId = "type",
                                      label = "Cancer type",
                                      choices = c("COAD","BRCA", "KIRC", "KIRP","BLCA","HNSC","KICH","LIHC","LUAD","LUSC","PRAD",
                                                  "STAD","THCA","UCEC")),
                          selectizeInput(inputId = "gene",
                                         label = "Lipid gene",
                                         choices =read.table(file = "lipidIDmatched.txt", header =TRUE),
                                         options = list(maxOptions = 3000)),
                          submitButton(text = "Submit")),
                        
                        
                        
                        
                        mainPanel(plotOutput("Survival"),
                                  downloadButton('downloadPlot', 'Download Plot')
                                  #tabPanel("Differential gene expression", plotOutput("expression"),value = "differential")
                                  #             tabPanel("vocanoPlot", plotOutput("vocanoPlot"))),
                                  
                                  
                        ))))),
    tabPanel("HELP",
             fluidRow(
               column(12, offset = 0.5,
                      img(src="TCGA_1.jpg", width=1900, height=644),
                      br(),
                      h4("In total, we have five pages deployed in the shiny interactive site. We show the data we rely on on the site and the methods we use for data output on page", strong("①"),"：Home. We introduce the use guide and related description of the website on page ", strong("③"),"HELP. On page ",strong("④"),":CITATION,We put the paper on this page, if you want to use the specific data and detailed data information of this website, you can click on this page to cite the paper",strong("⑤"),":ABOUTUS, we briefly introduce the people who contributed 
                    to the construction of this site and provide the source code and our contact information. If you have any questions or suggestions about this page, please feel free to contact us. The most important page in this website is page" , strong("②")," :Survival analysis page, where you can analyze the KM survival of 14 cancer species and their 
                    corresponding genetic data obtained by us. When you want to choose the type of cancer, you can change the type of cancer in button ", strong("⑥"),". When you want to change the gene, you can change the gene you want in button ", strong("⑦"),". Finally, you can click button ", strong("⑧"), "to confirm your choice of cancer and genes. Then you can see the KM survival
                    analysis curve corresponding to the cancer type and gene of your choice.Of course, if you want to download the survival analysis images of your selected cancer species and genes, you can click button  ",strong("⑨"),". so that you can download the images drawn on this website to your computer.")
               ))),
    tabPanel("CITATION",
             fluidRow(
               column(12, offset = 0.5, 
                      mainPanel(
                        br(),
                        br(),
                        h3("Please cite the following paper when using the server and data at this website:"),
                        h4("Guoqing Liu1,2#*, Yan Yang1#, Xuejia Kang1#, Hao Xu1, Jing Ai1, Min Cao1, Guojun Liu1,2*.A pan-cancer analysis of lipid metabolic alterations in primary and metastatic cancers.
                                                                                                                  (The paper is in the review stage for publication, and the data will be updated after publication)"),
                        
                        
                      )
                      
                      
                      
               ))),
    tabPanel("ABOUT US",
             fluidRow(
               column(12, offset = 0.5,
                      mainPanel(h3("About Us"),h4( "This is an interactive website developed by several bioinformatics professors and students from Inner Mongolia University of Science & Technology.The website is one of the exhibits for the team to participate in the National Life Science Competition for College Students (2022) and to publish related papers. 
                              Importantly,These results are wholly or partly based on the data results generated by TCGA research network. You can find the data we studied in",a(" TCGA database",href = "https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga"),strong("(https://www.cancer.gov/ccg/research/genome-sequencing/tcga)") ,"for analysis to prove the correctness of our results. "),
                                br(),
                                h3("The function:"), h4("At present, this website has 14 kinds of tumor sample information KM survival analysis image drawing and download function.If necessary in the later stage, we will improve and change it."),
                                br(),
                                h3("Questions:"), h4("The R source code for this site can be found on", a("Github", href="https://github.com/Whaoe-hao?tab=repositories"),strong("(https://github.com/Whaoe-hao?tab=repositories)"),",and if you still have questions, you can send an email below ."),
                                br(),
                                h3("Contact:"), h4("Please email with any questions, comments, or suggestions to", a("Guoqing Liu.", href="mailto:gqliu1010@163.com"),strong("(gqliu1010@163.com)"))
                      )
               )
             ))
  ))


server <- function(input,output){
  
  gene_id <- reactive({
    input$gene
    
  })
  
  cancer_data <- reactive({
    list(input$type,
         fileName <-paste0("R.DATA/",input$type,".RData"),
         load(fileName),
         datGroup <- get("datGroup"),
         my.surv <- get("my.surv")
         
         # 上面两行代码也可以用下面三行替代
         # load(fileName),
         # datGroup <- eval(parse(text = "datGroup")),
         # my.surv <- eval(parse(text = "my.surv"))
    )
  })
  
  
  
  output$Survival <- renderPlot({
    
    #############cancer type selection and read data###########################
    cancerData <- cancer_data()
    ##print(cancerData)
    datGroup <- cancerData[[4]]
    my.surv <-cancerData[[5]]
    
    
    
    #############KM survival analysis###########################
    geneName <- gene_id()
    # group=as.vector(t(datGroup[,gene_id()]))
    #kmfit1 <- survfit(my.surv~group,data=datGroup)
    
    
    kmfit1 <- survfit(as.formula(paste("my.surv~",geneName)),data=datGroup)
    
    
    # # 检验显著性
    diffResult=survdiff(as.formula(paste("my.surv~",geneName)),data=datGroup)
    pvalue=1-pchisq(diffResult$chisq,length(diffResult$n)-1)
    
    
    par(mai=c(2,3,0.5,0.5))
    par(cex=2)
    plot(kmfit1,las=1,lwd=3,col = c("coral","blue"),main=geneName,xlab = "Time (day)", ylab = "Survival probability")
    legend("topright", legend=c("High","Low"),col = c("coral","blue"), bty="n", lwd=3,xpd=T)
    text(x=10,y=0.1,labels=paste0("p=",round(pvalue,3)),adj=0)
    #ggsurvplot(kmfit1,data=datGroup,conf.int =F, pval = T,pval.coord=c(0.1,0.24),xlab="Time (day)",ylim=c(0.22,1),ggtheme=theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill='transparent',color='black'),legend.key =element_rect(fill=NA,color=NA),legend.background=element_rect(fill = 'NA',colour = 'NA'),legend.direction='vertical',legend.margin=margin(0, 0, 0, 0,"pt")),title="",legend.title = element_blank(),legend=c(0.85,0.88),legend.labs=c("High","Low")) #其他程序中运行此语句正常,但是在shiny中运行出错
    
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("survival_plot_", input$type, "_", input$gene, ".png", sep='')
    },
    content = function(file) {
      # 保存当前的图像
      png(file)
      cancerData <- cancer_data()
      datGroup <- cancerData[[4]]
      my.surv <-cancerData[[5]]
      geneName <- gene_id()
      kmfit1 <- survfit(as.formula(paste("my.surv~",geneName)),data=datGroup)
      diffResult=survdiff(as.formula(paste("my.surv~",geneName)),data=datGroup)
      pvalue=1-pchisq(diffResult$chisq,length(diffResult$n)-1)
      par(mai=c(2,3,0.5,0.5))
      par(cex=2)
      plot(kmfit1,las=1,lwd=3,col = c("coral","blue"),main=geneName,xlab = "Time (day)", ylab = "Survival probability")
      legend("topright", legend=c("High","Low"),col = c("coral","blue"), bty="n", lwd=3,xpd=T)
      text(x=10,y=0.1,labels=paste0("p=",round(pvalue,3)),adj=0)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)

