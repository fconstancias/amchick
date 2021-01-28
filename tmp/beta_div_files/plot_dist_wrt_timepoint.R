#' @title Plot beta diversity metrics relative to certain timepoint over time
#' @author Sneha Sundar and Florentin Constancias
#' @param 
#' 
#' @return 

#' @export
#' @examples 


phyloseq_plot_beta_div_wrt_timepoint <- function(distances, 
                                                           bdiv_list,
                                                           physeq,
                                                           timepoint,
                                                          group_var,
                                                          time_var,
                                                          group_to_compare=NULL
                                                           ){
  require(phyloseq)
  require(microbiome)
  require(tidyverse)
  require(usedist)
  
  #function name wrt timepoint
  
  #re
  
  
  plot_beta_div_wrt_timepoint <- function(dist, 
                                                    bdiv_list,
                                                    physeq,
                                                    timepoint, 
                                                    group_var,
                                                    time_var,
                                                    group_to_compare=NULL){
  
  
  
  #time ref= previous
  d.mat <- bdiv_list[[dist]]
  
  as(sample_data(physeq),"matrix") %>%
    data.frame(check.names=FALSE) %>%
    rownames_to_column("Sample_ID") -> sample.data
  
  sample.data[,time_var] <- as.numeric(sample.data[,time_var])
 
  
  stopifnot(all.equal(labels(d.mat), sample.data$Sample_ID))
  
  item_groups <- sample.data[,group_var]
  
  dist_df <- usedist::dist_groups(d.mat,item_groups)
  
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
    
    days.df <- days.df %>% 
               ungroup() %>% 
               mutate(Day2=day2)
    
    days.df.reverse <- days.df %>% 
                      rename(Day2=Day1,Day1=Day2)
    
    #combine both datarames to get the complete one
    days.df.complete <- rbind(days.df,days.df.reverse)
    
    days.df.complete[,group_var] <- paste("Within",pull(days.df.complete,.data[[group_var]]))
    
    #Getting the right dataframe for plotting
    df_plot <- semi_join(dist_df,days.df.complete,
                       by=c("Label"= group_var,"varGroup1"="Day1","varGroup2"="Day2")) %>% 
              arrange(varGroup1) %>% 
              arrange(varGroup2)
    
    df_plot %>% 
      group_by(Label) %>% 
      arrange(Label) %>% 
      arrange(varGroup1,.by_group=TRUE) %>% 
      arrange(varGroup2,.by_group=TRUE) -> df_plot
    
    #if day1 > day2 swap them 
    swap_indices<-which(df_plot$varGroup1 > df_plot$varGroup2)
    
    df_plot[swap_indices,c("varGroup1","varGroup2")] <- df_plot[swap_indices,c("varGroup2","varGroup1")]
    
    #plotting
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
    
    if(is.null(group_to_compare)){
      stop("Error: You need to specify the name of the common group to compare. This group needs to be one of the categories in your group_var argument.")
    }
    
    #keeping only between sample distances with the common group
    
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
    
    
    df_plot<-dist_df %>% group_by(Group2) %>% arrange(varGroup2,.by_group=TRUE)
    
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
      theme_bw() + xlab("Day") + ylab(paste("Distance to",group_to_compare)) -> plot
    
    
    
    
  }
  
  
  return(plot)
  
  }
  #you will use lapply to apply above function to a numver of distances
  
  
  
  res<- lapply(X=distances,FUN=plot_beta_div_wrt_timepoint,bdiv_list,
                physeq,
                timepoint, 
                group_var,
                time_var,group_to_compare)
  
  names(res) <- distances
  
  return(res)
  
}


