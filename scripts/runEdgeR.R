#!/project/chadm/Scott/Software/R/R-3.3.2/bin/Rscript --vanilla
#!/usr/bin/Rscript --vanilla
# The user may need to change their Rscript directory, but make sure
# that the --vanilla option is enabled to ensure reproducibility!
# Author: Sean McIlwain
# Edited by Scott Simpkins

# Installs edgeR if not currently installed
if (!require(edgeR)) {
    source('http://bioconductor.org/biocLite.R')
    biocLite('edgeR')
}

# Installs corrplot if not currently installed
if (!require(corrplot)) {
    r = getOption('repos')
    r['CRAN'] = 'http://cran.rstudio.com/'
    options(repos = r)
    install.packages('corrplot')
}



########################################
# FUNCTIONS AND PROCEDURES
########################################
usage<-function() {

  usage_string = paste0("usage() : runEdgeR.R <data file> <control file>\n",
  	       	        "  --threshold n : threshold to use for genes");

  return(usage_string);	  
}


doEdgeR<-function(data_in, col1, col2, col3=NA, threshold=0, file_stem=NA) {
  #source("http://www.bioconductor.org/biocLite.R")
  #biocLite("edgeR")
  library(edgeR)

  cat("doFoldChangesEdgeR for:",file_stem,"\n");
  data1 = data_in;
#  data1 = filterLowCountReads(data_in, col1, col2, FALSE, threshold);
  rownames(data1) = data1$b.no;
  data1 = data1[,c(col1,col2)];



  group = c(rep("Condition1", length(col1)), rep("Condition2", length(col2)));
  cds = DGEList(data1, group=group);
  cds = calcNormFactors(cds);  
  cds = estimateCommonDisp(cds);
  cds = estimateTagwiseDisp( cds );
  et = exactTest(cds);
  
  y = topTags(et, n=nrow(data1))$table;
  #print(y);
  ans = data.frame(
    b.no = rownames(y),
    foldChange=2^y$logFC,
    pvalue=y$PValue,
    padj = y$FDR,
    stringsAsFactors=FALSE);
  ans$padj[is.na(ans$padj)] = 1.0;
  ans$padj[is.na(ans$pvalue)] = 1.0;
  #ans = ans[!is.infinite(ans$foldChange) & !is.na(ans$padj),];

  return(ans);
}



getNumericColumnIndices<-function(data_in) {

  ans = c();
  for (col in 1:ncol(data_in)) {
    if (is.numeric(data_in[,col])) {
      ans = c(ans, col);
    }
  }
  return(ans);
}

getNumericColumns<-function(data_in) {
  return(data_in[,getNumericColumnIndices(data_in)]);
}


makeVolcanoPlots<-function(data_, fold_change_column, pvalue_column, title="", filestem=NA) {
  makeVolcanoPlot(data_, fold_change_column, pvalue_column, title);
  
  if (!is.na(filestem)) {
  pdf(makeExtFilePath(filestem,"volcano","pdf"))
  makeVolcanoPlot(data_, fold_change_column, pvalue_column, title);
  dev.off();
  
  png(makeExtFilePath(filestem,"volcano","png"));
  makeVolcanoPlot(data_, fold_change_column, pvalue_column, title);
  dev.off();
  }
  
  
  
}


makeVolcanoPlot<-function(data_, fold_change_column, pvalue_column, title_ = "Volcano!") {

ptsize=2;


pdata = data.frame(b.no=data_$b.no, x = data_[,fold_change_column], y = -log10(data_[,pvalue_column]));
pdata_ninf = pdata$y[!is.infinite(pdata$y)];
ninf = max(pdata_ninf);
pdata$y[is.infinite(pdata$y)] = ninf;

pdata_ninf = pdata$x[!is.infinite(pdata$x)];
ninf = max(pdata_ninf);
pdata$x[is.infinite(pdata$x)] = ninf;


#pdata = getEcoliBnosToNames(pdata);

#Color all points grey

x = max(max(pdata$x),1/min(pdata$x));

xlim_=c(1/x,x);
ylim_=c(0,max(4,max(pdata$y)));

print(xlim_);
print(ylim_);

plot(x=pdata$x,
    y=pdata$y,
    main=title_,
    xlab="Fold-Change",
    ylab="-log10(p-value)",
    log="x",
    xlim=xlim_,
    ylim=ylim_,
    col="grey",pch=19,lwd=max(1,ptsize-2));

#Color all points with pvalue < 0.05 dark grey
sig_data = pdata[pdata$y>-log10(0.05),];
points(x=sig_data$x,y=sig_data$y,pch=19,col="darkgrey", lwd=ptsize);

lines(xlim_, c(-log10(0.05), -log10(0.05)),col="blue"); # 5% FDR
lines(c(0.5,0.5), ylim_, col = "red"); #1/2 fold change
lines(c(2, 2), ylim_, col="red"); # 2 -fold change

#makeVolcanoPlotC(pdata, gene_colors$form_detox_genes, gene_colors$form_detox_color, ptsize)
#makeVolcanoPlotC(pdata, "rpoS", rpos_color, ptsize);


  legend(2,ylim_[1]+1, pch = 19, cex=0.841,
    legend=c("adj. p-value < 0.05"),
    col=c("darkgrey"));

  return(sig_data);

}


############################
# makeFilePath - takes the date_tag, a globably define variable 
# the filestem, and the end to build a complete path
############################
makeFilePath<-function(prefix, file, postfix) {

  if (prefix != "") {
    if (file != "") {
      ans=paste0(prefix,".",file);
    } else {
      ans = prefix;
    }
  } else if (file != ""){
    ans=file;
  } else {
    error("both file and prefix are blank!");
  }
  
  if (exists("date_tag")) {
    ans = paste0(ans, ".", date_tag);
  } else {
    #cat("Date tag doesn't exist!\n");
  }
  ans = paste0(ans, ".", postfix);
  return(ans);

}



###########################
# This function will create a file path with a
# subdirectory in ext and of the format ext/prefix.file_stem_.ext

makeExtFilePath<-function(prefix,file_stem_,ext_) {

  dir.create(ext_, showWarnings=FALSE);

  if (prefix == "") {
    ans = makeFilePath("",paste0(ext_,"/",file_stem_), ext_);
  } else {
    ans=makeFilePath(paste0(ext_,"/",prefix),file_stem_,ext_);
  }
  cat(ans,"\n");
  return(ans);
}



filterLowCountReadsMax<-function(data_in, threshold) {

    data_num = getNumericColumns(data_in);
    max_count = apply(data_num,1,max, na.rm=TRUE)
    
    #hist(max_count, breaks=max(max_count))
    
    nfilter = sum(max_count > threshold)
    
    cat("Keeping ",nfilter," out of ",nrow(data_in),"\n");
    ans = data_in[max_count > threshold,];
    
    return(ans);

}

###############
# Reads in the data file
###############
loadData<-function(path, tag=NA) {

  ptr = file(path);
  header_line = readLines(ptr, 1);
  close(ptr);
  header_tokens = strsplit(header_line,"\t")[[1]];
  #print(header_tokens);
  ans = read.csv(path, stringsAsFactors=FALSE, sep="\t");
  colnames(ans) = header_tokens;
  print(colnames(ans));
  colnames(ans)[1] = 'ORF'
  ans$ORF = vapply(strsplit(ans[['ORF']], '_'), `[`, character(1), 1)
  
  #ans = ans[ans$ORF!="EWEIGHT",];
  #ans = ans[,colnames(ans) != "GWEIGHT"];
  #ans = ans[,colnames(ans) != "NAME"];
  
  if (!is.na(tag)) {
    colnames(ans)[2:ncol(ans)] = paste0(tag,".",colnames(ans)[2:ncol(ans)]);
  }
  return(ans);    

}

###################
# Reads in the list of controls from the control file
###################
loadControls<-function(path) {
  
  ptr = file(path);
  controls = readLines(ptr);
  close(ptr);
  for (idx in 1:length(controls)) {
    if (isDigit(controls[idx],1)) {
        warning("control ",controls[idx], " starts with a digit!\n");
        #stop();
    }
  }
  
  return(controls)    
}


getBnosToOrfs<-function(data_in) {
    ans = rep(NA, nrow(data_in));
    for (idx in 1:nrow(data_in)) {
        bno = data_in$b.no[idx];
        orf = bno_to_orf_df$ORF[bno_to_orf_df$b.no == bno];
        ans[idx] = orf;
    }
    return(ans);
}


makeBoxPlots<-function(data_, filestem_ = NA, title_ = "", ylabel = "log2(Intensity)") {

  if (!is.na(filestem_)) {

    png(makeExtFilePath(filestem_,"box","png"));
    makeBoxPlot(data_, title_, ylabel);
    dev.off();
    pdf(makeExtFilePath(filestem_,"box","pdf"));
    makeBoxPlot(data_, title_, ylabel);
    dev.off();
  } else {
    makeBoxPlot(data_, title_, ylabel);
  }

}

makeBoxPlot<-function(data_, title_ = "", ylabel = "log2(Intensity)") {
  par(las=2);
  par(mar=c(10,4,4,2))
  boxplot(getNumericColumns(data_), ylab="log2(Intensity)", title=title_, cex.axis=0.50);


}
########################################
#Checks whether the character in the indexth part of the string is a 
#number
########################################
isDigit<-function(in_string, index=1) {

  s = substr(in_string, index, index);
  
  ans = s == "0";
  ans = ans | s == "1";
  ans = ans | s == "2";
  ans = ans | s == "3";
  ans = ans | s == "4";
  ans = ans | s == "5";
  ans = ans | s == "6";
  ans = ans | s == "7";
  ans = ans | s == "8";
  ans = ans | s == "9";
  
  return(ans);

}


log2Data<-function(data_, shift_=0) {

  ans = data_;

  for (col_idx in 1:ncol(data_)) {
    if (is.numeric(data_[,col_idx])) {
      ans[,col_idx] = log2(data_[,col_idx]+shift_);
    } else {
      ans[,col_idx] = data_[,col_idx];
    }
  }
  return(ans);
}


##################################
# CreateCorrelationPlot
##################################
createCorrelationPlots<-function(data, filestem, method="pearson") {

  png(file=makeExtFilePath(filestem,"corr","png"));
  createCorrelationPlot(data, method);
  dev.off();

  pdf(file=makeExtFilePath(filestem,"corr","pdf"));
  createCorrelationPlot(data, method);
  dev.off();

  print(str(getNumericColumns(data)))
  correlation=cor(getNumericColumns(data), use="pairwise.complete.obs");
  print(str(correlation))
  write.csv(correlation, makeExtFilePath(filestem, "corr", "csv"));		
  return(correlation);	 

}


createCorrelationPlot<-function(data, method="pearson") {
  col1 = colorRampPalette(c('black','black','black','red','red','orange','yellow','green'))(200);
  correlation=cor(getNumericColumns(data), use="pairwise.complete.obs", method=method);
  NAs = is.na(correlation);
  if (sum(as.vector(NAs)) > 1) {
    cat("Warning ", sum(NAs), " detected in correlation\n");
    correlation[NAs] = 0.0;

  }
  
  library(corrplot)
  corrplot(
    correlation, 
    method="color", 
    tl.cex=0.4,
    cl.cex=0.5,
    col=col1,
    addCoef.col='black'
    #cl.lim=c(0,1)
);

}


##############################
# MAIN SCRIPT
##############################

args = commandArgs(trailingOnly = TRUE);

data_path = NA;
control_path = NA;
threshold = -1;

arg_idx = 1;
main_args = c();
while (arg_idx <= length(args)) {
  if (args[arg_idx] == "--threshold") {
     arg_idx = arg_idx + 1;
     threshold = as.integer(args[arg_idx]);
  } else {
    main_args = c(main_args, args[arg_idx]);
  }
  arg_idx = arg_idx + 1;
}

if (length(main_args) != 2) {
   stop(usage())
}

data_path = main_args[1];
control_path = main_args[2];


cat("Threshold:", threshold, "\n");
cat("Data:", data_path,"\n");
cat("Control:",control_path,"\n");


#Load data

raw_data = loadData(data_path);

controls = loadControls(control_path);

print(str(raw_data))
print(str(controls))

cat("=======================\n");
cat("CONTROLS\n");
print(controls);
cat("=======================\n");

bno_to_orf_df = data.frame(b.no = 1:nrow(raw_data), ORF=raw_data$ORF, stringsAsFactors=FALSE);

#Sort by controls first

conditions = rep(NA, ncol(raw_data));
for (idx in 2:ncol(raw_data)) {
    temp = strsplit(colnames(raw_data)[idx], "_")[[1]];
    condition = paste(temp[1:(length(temp)-1)], sep="_",collapse="_");
    conditions[idx] = condition;
}

uconditions = table(na.omit(conditions));

cat("Condition Table:\n");
print(uconditions);

data_new = data.frame(b.no=as.character(1:nrow(raw_data)));
data_new = cbind(data_new, raw_data[,(colnames(raw_data) %in% controls)]);

conditions_new = c(NA, rep("Control", length(controls)));


for (idx1 in 1:length(uconditions)) {
    condition = names(uconditions)[idx1];
    for (idx2 in 2:ncol(raw_data)) {

    	current_name = colnames(raw_data)[idx2];
	current_condition = conditions[idx2];

	if ((!current_name %in% controls) & current_condition == condition) {
	   data_new = cbind(data_new, raw_data[idx2]);
	   conditions_new = c(conditions_new, condition);
	}
    }
}

conditions = conditions_new;
uconditions = unique(na.omit(conditions));

data_filtered = filterLowCountReadsMax(data_new, threshold);

if (min(getNumericColumns(data_filtered)) == 0) {

   data_filtered[,getNumericColumnIndices(data_filtered)] = 
     data_filtered[,getNumericColumnIndices(data_filtered)] + 1;

}


rownames(data_filtered) = data_filtered$b.no;
write.csv(data_filtered, makeExtFilePath("raw", "filtered", "csv"), row.names=FALSE);

library("edgeR")

print(colnames(data_filtered));
cat("========================\n");
print(conditions);

data_cds = DGEList(data_filtered[,!is.na(conditions)], group = conditions[!is.na(conditions)]);
data_cds = calcNormFactors(data_cds);
data_cpm = data.frame(cpm(data_cds));
data_cpm$b.no = data_filtered$b.no;
data_cpm$ORF = getBnosToOrfs(data_cpm);

makeBoxPlots(log2Data(data_filtered, 1), "raw", "Raw Counts");
makeBoxPlots(log2Data(data_cpm, 1), "cpm", "Counts Per Million");

pdf(makeExtFilePath("cpm","pairs","pdf"));

try(pairs(getNumericColumns(data_cpm)));
dev.off();

str(data_filtered)
str(data_cpm)

corr_raw = createCorrelationPlots(data_filtered, "raw");
corr_cpm = createCorrelationPlots(data_cpm, "cpm");
write.csv(corr_cpm, makeExtFilePath("corr","cpm","csv"), row.names=FALSE);

############
# Pairwise Comparisons
############

cond1_c = controls;

data_fc = data.frame(b.no = data_filtered$b.no);    
for (idx1 in 1:length(uconditions)) {
  
  if (!is.na(uconditions[idx1]) && uconditions[idx1] != "Control") {
     cond2_c = colnames(data_filtered)[!is.na(conditions) & conditions == uconditions[idx1]];
     edger_fc = doEdgeR(data_filtered, cond1_c, cond2_c, threshold=-1);
     makeVolcanoPlots(
       edger_fc, 
       "foldChange", 
       "padj", 
       paste0(uconditions[idx1],".v.","Control"), paste0(uconditions[idx1],".v.","Control"));

     data_fc = merge(data_fc, edger_fc[,c("b.no","foldChange")], all.x=TRUE, all.y=TRUE, by.x="b.no", by.y="b.no");
     colnames(data_fc)[ncol(data_fc)] = uconditions[idx1];
     orfs = getBnosToOrfs(edger_fc);
     edger_fc = cbind(orfs, edger_fc);
     colnames(edger_fc)[1] = "ORF";
     edger_fc = edger_fc[,c("ORF","foldChange","pvalue","padj")];
     write.csv(edger_fc, makeExtFilePath(paste0(uconditions[idx1],".v.","Control"), "pw.rc.edger","csv"), row.names=FALSE);
     gene_lists = as.character(edger_fc$ORF[edger_fc$padj < 0.05]);
     cat("Selected ",length(gene_lists)," out of ",nrow(edger_fc),"\n");
  }
  cat(idx1, " out of ", length(uconditions), " conditions\n");
}











