# 如果需要的library没有的话，先下载
if (!"optparse" %in% installed.packages()[, "Package"]) {
  # 这个包用来解析命令行选项
  install.packages("optparse", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
library(optparse)
# 解析第一个选项
args <- commandArgs(trailingOnly = TRUE)
command <- args[1]
library("htmltools")
library("webshot")    
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2){
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}


get_klebsiella_species <- function(contigs, datafolder) {
  # 执行mash dist命令
  system("mkdir temp")
  mash_command <- paste("mash dist", file.path(datafolder, "species_mash_sketches.msh"), contigs)
  mash_output <- system(mash_command, intern=TRUE)
  mash_command2 <- paste0("mash sketch -p 2 -o ",contigs," ",contigs)
  system(mash_command2, intern=TRUE)
  best_species <- NULL
  best_distance <- 1.0
  
  # 处理命令的输出
  for (line in mash_output) {
    line_parts <- strsplit(line, "\t")[[1]]
    reference <- line_parts[1]
    if (length(line_parts) < 4) {
      next
    }
    species <- strsplit(reference, "/")[[1]][1]
    distance <- as.numeric(line_parts[3])
    
    # 处理物种名称
    species <- gsub("Raoultella", "Klebsiella (Raoultella)", species)
    species <- gsub("_", " ", species)
    species <- gsub(" subsp ", " subsp. ", species)
    if (grepl(" unknown$", species)) {
      species <- gsub(" unknown", " (unknown species)", species)
    }
    
    if (distance < best_distance) {
      best_distance <- distance
      best_species <- species
    }
  }
  
  # 返回最佳物种和置信度
  if (best_distance <= 0.02) {
    return(list(best_species, "strong"))
  } else if (best_distance <= 0.04) {
    return(list(best_species, "weak"))
  } else {
    return(list("unknown", ""))
  }
}
if (command == "plot") {
  suppressMessages(library(ggtree))
  library(phytools)
  library(formattable)
  suppressMessages(library(tidyverse))
  library(jsonlite)  
  library(readr)
  library(reshape2)
  #命令行选项
  option_list <- list(
    make_option(c("--config"), type = "character", default = NULL, help = "Path to JSON config file", metavar = "character"),
    make_option(c("--output"), type = "character", default = "output.png", help = "Output file name with extension", metavar = "character")
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(opt_parser, args = args[-1])
  json_file_path <- opts$config
  data_file_path <- opts$input
  output_file_name <- opts$output
  create_heatmp_plot_from_json <- function(json_file_path, data_file_path, output_file_name) {
    #####文件夹切换#####
    if (!file.exists(json_file_path)) {
      stop(paste("JSON file not found:", json_file_path))
    }
    
    tryCatch(
      {
        config <- fromJSON(json_file_path)
      },
      error = function(e) {
        stop(paste("Error reading JSON file:", e))
      }
    )
    # 从JSON文件中读取列配置
    column_config <- config$columns
    
    # 从JSON中读取数据文件path
    data_file_path <- config$dataFile$dataFilePath
    
    # 检查CSV文件是否存在
    if (!file.exists(data_file_path)) {
      stop(paste("Path file not found:", data_file_path))
    }
    
    tryCatch(
      {
        data <- read.csv(data_file_path)
      },
      error = function(e) {
        stop(paste("Error reading Path file:", e))
      }
    )
    
    
    ######检查软件
    ######检查是否存在Path列
    if (!column_config$path %in% colnames(data)) {
      stop("path col not found in data")
    }
    colnames(data)[colnames(data) == column_config$path] <- "path"
    
    
    ######把文件挪到一个文件夹下
    
    system("mkdir -p temp")
    for (i in row.names(data)) {
      tryCatch({
        temp <- paste0("cp ", data[i, "path"], "/*.fasta temp/")
        system(temp)
      }, error = function(e) {
        # 处理错误的操作
        print(paste("Error occurred for row", data[i, "path"]))
        # 抛出错误并终止程序
        stop("Error occurred. Program terminated.")
      })
    }
    ####检查Fasta数量，小于3 绘制表格 大于3绘制进化树
    fasta_shuliang=list.files("./temp",pattern = "*fasta$")
    
    if (length(fasta_shuliang) >2) {
      print("File > 3 ,print tree plot")
    } else {
      print("File < 3  ,print table")
    }
    if (length(fasta_shuliang) >2) {
    mode="tree"
    } else {
     mode="table"
    }
    if (mode=="tree"){if (column_config$mode=="table") mode="table"}##fasta数量大于2，想要表格形式
    data_result=data.frame(label=c(),species=c(),strength=c())
    n=1
    for (fasta_file in fasta_shuliang) {
      fasta_file1=paste0("./temp/",fasta_file)
      data_temp=get_klebsiella_species(fasta_file1,column_config$databasepath)
      data_result[n,1]=fasta_file
      data_result[n,2]=data_temp[1]
      data_result[n,3]=data_temp[2]
      n=n+1
    }
    colnames(data_result)=c("label","species","strength")
    if (mode=="tree"){
      system("mash paste combined_sketch.msh  temp/*.msh")##生成总mash文件
      system("mash dist combined_sketch.msh combined_sketch.msh > combined_distance.tsv")#计算距离
      disfile<-"combined_distance.tsv"
      datdis  <- read.table(disfile,header=FALSE, stringsAsFactors =FALSE)
      datdis
      datdis$V1=gsub("./temp/","",datdis$V1)#去除temp
      datdis$V2=gsub("./temp/","",datdis$V2)
      colnames(datdis)  <- c("anim1","anim2","distr","comp4","comp5")
      
      datdis %>% 
        mutate(anim1c=(anim1),
               anim2c=(anim2)) -> datdis
      
      
      datsel  <- datdis  %>% select(anim1c,anim2c, distr)
      datwide  <- datsel  %>% pivot_wider(names_from = anim2c, values_from = distr)
      
      datmat  <- as.matrix(datwide  %>% select(-anim1c))
      rownames(datmat)  <- datwide$anim1c
      
      datmat
      
      tr  <- nj(datmat)
      
      #group_file
      groupInfo_CAZR <- split(data_result$label, data_result$species)
      #重要，必须两步midpoint
      tr <- groupOTU(tr, groupInfo_CAZR)
      groupInfo_CAZR <- split(data_result$label, data_result$strength)
      tr <- groupOTU(tr, groupInfo_CAZR,group_name="strength")
      
      test=ggtree(tr,
                  layout=config$plot_settings$layout,
                  branch.length=config$plot_settings$branchlength,
                  size=config$plot_settings$treesize,
                  color=config$plot_settings$treecolor,
                  linetype=config$plot_settings$treelinetype)
      
      ####title设定###############
      title_temp=labs(title=config$general$title)
      
      ########标尺设定############
      scale1=theme()
      if (config$general$treescaletype=="xasis")scale1= theme_tree()
      fontsize=config$plot_settings$scalefontsize
      linesize=config$plot_settings$scalelinesize
      offset=config$plot_settings$scaleoffset
      if (config$general$treescaletype=="dandu")scale1=geom_treescale(fontsize=fontsize,
                                                                      linesize=linesize,
                                                                      offset=offset)
      
      
      ########tiplab设定(theme()无任何效果)#########
      if(config$general$tiplab){
        temp_tip=geom_tiplab(size=config$plot_settings$tiplabsize)
        temp_tip2=geom_tiplab(aes(label=group),size=config$plot_settings$tiplabsize,offset=config$plot_settings$staroffset)
        temp_tip3=geom_tiplab(aes(label=strength),size=config$plot_settings$tiplabsize,offset=config$plot_settings$baroffset)
      }else{temp_tip=theme()
      temp_tip2=theme()
      temp_tip3=theme()
        }
      test+title_temp+temp_tip+temp_tip2+temp_tip3+scale1
      ggsave(output_file_name,height = config$plot_settings$height,width = config$plot_settings$width)
    }
    if (mode=="table"){
      FT <- formattable(data_result)
      
      export_formattable(FT,output_file_name)
      
     
    }
    
  
  }
  create_heatmp_plot_from_json(json_file_path,data_file_path, output_file_name)
}
