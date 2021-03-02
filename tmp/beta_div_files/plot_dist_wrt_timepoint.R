#' @title Plot time evolution of beta diversity metrics with respect to specified reference timepoints 
#' @author Sneha Sundar and Florentin Constancias
#' @param distances a character vector of the names of the distance metrics to be plotted; valid distances should be objects in \code{bdiv_list}
#' @param bdiv_list the output of function \code{phyloseq_compute_bdiv} containing a list of computed distance matrices
#' @param physeq rarefied phyloseq object
#' @param timepoint specifies which distances should be plotted over time; can be "previous", "fixed" or "between.ref.group". More details on what these mean in notes.
#' @param group_var name of the metadata column containing group information
#' @param time_var name of the metadata column containing time information
#' @param group_to_compare name of the reference group for comparison of distances. Should be specified only if \code{timepoint} is "between.ref.group" 
#' @param fixed_time the fixed time to compare distances within each group. Should be specified only if \code{timepoint} is "fixed" 
#' @return a list of ggplot objects named according to the distance metric that was plotted; the plots show us how a specified kind of beta diversity metric evolves over time for each group.
#' @note 
#' If \code{timepoint} is "previous",within each group in \code{group_var} , 
#' we pick out all the distances comparing current timepoint to previous 
#' timepoint and plot them over time. 
#' 
#' If \code{timepoint} is "between.ref.group", we pick out all the distances 
#' between group_to_compare and group in code{group_var} except \code{group_to_compare} 
#' for which the times are the same. A simple example: if "A", "B","C" are the groups
#' and "A" is set as \code{group_to_compare}, then distances between Day.1 of B and Day.1 of A, 
#' Day.2 of B and Day.2 of A,.. etc  are plotted for each group (we get one panel for "B" and another for "C"). 
#' 
#' If \code{timepoint} is "fixed", then within each group , distances compared 
#' to the \code{fixed_time} are plotted over time. 
#' 
#' 
#' 
#' @examples 
#' phyloseq_plot_beta_div_wrt_timepoint(distances = c("bray","bjaccard","wjaccard"),
#'                                    bdiv_list,
#'                                    physeq=ps_polyFermS_rare, 
#'                                    timepoint="between.ref.group",
#'                                    group_vartime_var="Day_from_Inoculum",                                        
#'                                    group_to_compare="CR_UNTREATED")
#' 
#' 


phyloseq_plot_beta_div_wrt_timepoint <- function(distances, 
                                                 bdiv_list,
                                                 physeq,
                                                 timepoint,
                                                 group_var,
                                                 time_var,
                                                 group_to_compare=NULL,
                                                 fixed_time=NULL
){
  require(phyloseq)
  require(microbiome)
  require(tidyverse)
  require(usedist)
  
  
  
  
  plot_beta_div_wrt_timepoint <- function(dist, 
                                          bdiv_list,
                                          physeq,
                                          timepoint, 
                                          group_var,
                                          time_var,
                                          group_to_compare=NULL){
    
    # PREPROCESSING 
    
    #extract distance matrix of class dist
    d.mat <- bdiv_list[[dist]]
    
    #sample data as dataframe
    as(sample_data(physeq),"matrix") %>%
      data.frame(check.names=FALSE) %>%
      rownames_to_column("Sample_ID") -> sample.data
    
    #make time column is numeric 
    sample.data[,time_var] <- as.numeric(sample.data[,time_var])
    
    #check if the labels match, if they don't throw an error
    stopifnot(all.equal(labels(d.mat), sample.data$Sample_ID))
    
    #extract metadata column containing group information
    item_groups <- sample.data[,group_var]
    
    #use library usedist to make a neat dataframe that shows for each distance in 
    #`d.mat` which samples were used for calculation and the group info of both samples
    dist_df <- usedist::dist_groups(d.mat,item_groups)
    
    #add the time info of the samples to `dist_df` . this is what we will use for 
    #picking out the distances we want to plot
    meta_df <- sample.data %>% select(Sample_ID,.data[[time_var]])
    
    left_join(dist_df,
              meta_df %>%
                dplyr::rename("varGroup1" = .data[[time_var]]
                ),
              by = c("Item1" = "Sample_ID")) %>%
      left_join(meta_df %>%
                  dplyr::rename("varGroup2" = .data[[time_var]]),
                by = c("Item2" = "Sample_ID")) -> dist_df
    
    
    if(timepoint=="previous")
    {
      #only within group distances are needed now
      dist_df %>%
        dplyr::filter(grepl("Within", Label)) -> dist_df
      
      
      #Create a dataframe specifying the days that we need to filter from the distance dataframe
      
      sample.data %>% 
        select(Sample_ID,.data[[group_var]],Day1=.data[[time_var]]) %>% 
        group_by(.data[[group_var]]) %>% 
        arrange(Day1,.by_group=TRUE) -> days.df
      
      #the last day of each group should not be compared the first day of next group
      na_fills<-cumsum(days.df %>% group_size())
      #specify the reverse order so all the right comparisons are picked out
      day2<-c(days.df$Day1[-1],NA)
      day2[na_fills] <- NA
      
      #ungrouping adding Day2 column
      days.df <- days.df %>% 
        ungroup() %>% 
        mutate(Day2=day2)
      
      #specifying the opposite too just in case `dist_df` has it in this order
      #remmember distance between A and B is same as distance between B and A
      days.df.reverse <- days.df %>% 
        rename(Day2=Day1,Day1=Day2)
      
      #combine both dataframes to get the complete one containing all possible day combinations of interest to us
      days.df.complete <- rbind(days.df,days.df.reverse)
      
      #making group info same as that of Label column in `dist_df`
      days.df.complete[,group_var] <- paste("Within",pull(days.df.complete,.data[[group_var]]))
      
      #Getting the right dataframe for plotting: pick out all the distances specified by `days.df.complete`
      # from dist_df
      df_plot <- semi_join(dist_df,days.df.complete,
                           by=c("Label"= group_var,"varGroup1"="Day1","varGroup2"="Day2")) %>% 
        arrange(varGroup1) %>% 
        arrange(varGroup2)
      
      #arrange and group by label
      df_plot %>% 
        group_by(Label) %>% 
        arrange(Label) %>% 
        arrange(varGroup1,.by_group=TRUE) %>% 
        arrange(varGroup2,.by_group=TRUE) -> df_plot
      
      #if day1 > day2 swap them to keep a consistent order: Day2 > Day1
      swap_indices<-which(df_plot$varGroup1 > df_plot$varGroup2)
      
      df_plot[swap_indices,c("varGroup1","varGroup2")] <- df_plot[swap_indices,c("varGroup2","varGroup1")]
      
      #plotting a connected scatterplot and faceting by Label
      df_plot %>%
        ggplot(aes(x = as.factor(varGroup2) , y = Distance)) +
        geom_point(size=1, alpha=0.6, aes(colour = Label, group=Label),
                   position=position_jitterdodge(dodge.width=0.9)) + 
        geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),
                  size = 0.4, linetype = "dashed", inherit.aes = TRUE, aes(group=Label, color = Label),
                  position=position_jitterdodge(dodge.width=0.9)) +
        # geom_boxplot(aes(group = varGroup2 %>% as.factor()),
        # fill = "transparent",
        # outlier.colour = NA,alpha=0.4) +
        facet_grid(Label ~ ., scales = "fixed") +
        # ggrepel::geom_text_repel(cex=2,
        #                      aes(label= Group1),
        #                      segment.color = 'black',
        #                      segment.size = 0.5,
        #                      # nudge_x =  -4,
        #                      # nudge_y = 0,
        #                      df_all_t0 %>% filter(Day2 == 6 & metric=="wunifrac") %>% unique()) +
        theme_bw() + xlab("Day") + ylab("Distance to previous timepoint") -> plot
      
    }
    
    
    if(timepoint=="between.ref.group"){
      
      #throw error if group to compare is not specified
      if(is.null(group_to_compare)){
        stop("Error: You need to specify the name of the common group to compare. This group needs to be one of the categories in your group_var argument.")
      }
      
      #keeping only between sample distances with the common group. The days to be compared should be the same
      
      dist_df %>% dplyr::filter(!grepl("Within", Label)) ->dist_df
      
      dist_df %>% 
        filter(Group1==group_to_compare | Group2==group_to_compare) %>% 
        filter(varGroup1==varGroup2) -> dist_df
      
      
      #to keep a consistent format. reference group is group1 and the other is group2
      
      if(sum(dist_df$Group2==group_to_compare)!=0){
        swap_indices<-which(dist_df$Group2==group_to_compare)
        dist_df[swap_indices,c("Item1","Item2","Group1","Group2")] <- dist_df[swap_indices,c("Item2","Item1","Group2","Group1")]
        dist_df$Label <- paste("Between",dist_df$Group1,"and",dist_df$Group2)
        
        
      }
      
      #the dataframe we will use for plotting
      df_plot<-dist_df %>% group_by(Group2) %>% arrange(varGroup2,.by_group=TRUE)
      
      #plotting a connected scatterplot and faceting by Label=
      df_plot %>%
        ggplot(aes(x = as.factor(varGroup2) , y = Distance)) +
        geom_point(size=1, alpha=0.6, aes(colour = Label, group=Label),
                   position=position_jitterdodge(dodge.width=0.9)) + 
        geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),
                  size = 0.4, linetype = "dashed", inherit.aes = TRUE, aes(group=Label, color = Label),
                  position=position_jitterdodge(dodge.width=0.9)) +
        # geom_boxplot(aes(group = varGroup2 %>% as.factor()),
        # fill = "transparent",
        # outlier.colour = NA,alpha=0.4) +
        facet_grid(Label ~ ., scales = "fixed") +
        # ggrepel::geom_text_repel(cex=2,
        #                      aes(label= Group1),
        #                      segment.color = 'black',
        #                      segment.size = 0.5,
        #                      # nudge_x =  -4,
        #                      # nudge_y = 0,
        #                      df_all_t0 %>% filter(Day2 == 6 & metric=="wunifrac") %>% unique()) +
        theme_bw() + xlab("Day") + ylab(paste("Distance to ",group_to_compare)) -> plot
      
      
      
      
    }
    
    
    
    if(timepoint == "fixed"){
      
      #throw error if fixed_time is not specified
      if(is.null(fixed_time)){
        stop("Error: You need to specify the fixed time within each group using which we will pick out distances. This time needs to be one of the times in the `time_var` column that is common to all groups.")
      }
      
      
      
      #only within group distances are needed now
      dist_df %>%
        dplyr::filter(grepl("Within", Label)) -> dist_df
      
      #get all the distances computed with a sample of `fixed_time`
      dist_df %>% 
        filter(varGroup1==fixed_time | varGroup2==fixed_time) -> dist_df
      
    
      #to keep a consistent format: varGroup1 is the fixed time
      if(sum(dist_df$varGroup2==fixed_time)!=0){
        swap_indices<-which(dist_df$varGroup2==fixed_time)
        dist_df[swap_indices,c("Item1","Item2","Group1","Group2","varGroup1","varGroup2")] <-   dist_df[swap_indices,c("Item2","Item1","Group2","Group1","varGroup2","varGroup1")]
        
      }
      
      #group by label and arrange by day2
      dist_df %>% group_by(Label) %>% arrange(varGroup2,.by_group=TRUE) -> df_plot
      
      
      #plot
      df_plot %>%
        ggplot(aes(x = as.factor(varGroup2) , y = Distance)) +
        geom_point(size=1, alpha=0.6, aes(colour = Label, group=Label),
                   position=position_jitterdodge(dodge.width=0.9)) + 
        geom_path(arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),
                  size = 0.4, linetype = "dashed", inherit.aes = TRUE, aes(group=Label, color = Label),
                  position=position_jitterdodge(dodge.width=0.9)) +
        # geom_boxplot(aes(group = varGroup2 %>% as.factor()),
        # fill = "transparent",
        # outlier.colour = NA,alpha=0.4) +
        facet_grid(Label ~ ., scales = "fixed") +
        # ggrepel::geom_text_repel(cex=2,
        #                      aes(label= Group1),
        #                      segment.color = 'black',
        #                      segment.size = 0.5,
        #                      # nudge_x =  -4,
        #                      # nudge_y = 0,
        #                      df_all_t0 %>% filter(Day2 == 6 & metric=="wunifrac") %>% unique()) +
        theme_bw() + xlab("Day") + ylab(paste("Distance to",time_var,fixed_time)) -> plot
      
      
    }
    
    return(plot)
    
  }
  
  
 
 
  #using lapply to generate plots for a vector of distances
  res<- lapply(X=distances,FUN=plot_beta_div_wrt_timepoint,bdiv_list,
               physeq,
               timepoint, 
               group_var,
               time_var,group_to_compare)
  
  #names of the plots will be the name of the distance metric used to compute the distance matrix
  names(res) <- distances
  
  return(res)
  
}