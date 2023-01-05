#The script was simplified from the TeloBase version given it does not contain interactivity necessary in the web application (inputs replaced with a static choice and conditional branching whether to show modal window was removed,
#parts of the script for input checking remained in and are commented upon within the code)
library(metacoder) #IMPORTANT - requires 0.3.5 version due to reintegration of taxa functionality (available at https://github.com/grunwaldlab/metacoder/releases)
library(stringi)
library(plotly)

#save input (in TeloBase script this depends on the input from the taxonomy filter options)
kingdom<- NULL
phylum<- NULL
class<- NULL
order<- "Brassicales"
family<- NULL
genus<- NULL
species<- NULL

#create vector with vector lengths to check if no filter is filled and later which lowest taxonomy rank is actually filled
taxonomy_rank_length<- c(length(species), length(genus), length(family), length(order), length(class), length(phylum), length(kingdom))

#position of taxonomy_rank_length that has non-null value to find filter with lowest taxonomy rank that is actually filled
which_taxonomy_rank<- which(taxonomy_rank_length > 0)[1]

#load datatable based on the filtered taxonomy option (in TeloBase it is dependent on the taxonomy filter options, here statically uploads telomere_sequences_Brassicales.csv)
datafile_raw <- read.csv2("telomere_sequences_Brassicales.csv", stringsAsFactors = FALSE) 

#create taxonomy_rank vector to match which_taxonomy_rank to find the preciese taxonomy rank that was taken
taxonomy_rank<- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

#set which taxonomy rank was taken last and  match it to the taxonomy_rank
lowest_selected_rank<- taxonomy_rank[which_taxonomy_rank]

#create list from the inputs of taxonomy filters to use which_taxonomy_rank to get the precise vector
taxonomy_rank_list<- list(species, genus, family, order, class, phylum, kingdom)

#last selected taxonomy rank using which_taxonomy_rank
lowest_rank_name<- taxonomy_rank_list[[which_taxonomy_rank]]

#vector for checking if there is more than 1 selected element in the rank
length_lowest_rank_name<- length(lowest_rank_name)

#checker if more elements in the lowest rank have the same taxonomy (if not later switched to FALSE, it will give a modal back)
taxonomy_checker<- TRUE

datafile_highlight<- data.frame(to_highlight = 0, datafile_raw) #prepare setup for the sequence highlight

metacoder_sequence_to_plot<- "TTTAGGG" #(in case of TeloBase it is based on the input)

species_to_highlight<- datafile_raw$species[datafile_raw$sequence %in% metacoder_sequence_to_plot] #select species that do have metacoder_sequence_to_plot
datafile_highlight$to_highlight[datafile_highlight$species %in% species_to_highlight]<- 1 #change to_highlight for these species to 1 value

datafile_highlight<- datafile_highlight[,-c(3:8)] #remove uninportant columns to speed up the calculation
unique_df_species<- datafile_highlight[!duplicated(datafile_highlight$species),] #lower the amount of information for later plotting

#code in if changes lowest_selected_rank and lowest_rank_name accordingly
if (length_lowest_rank_name > 1) {
  
  #check if two kingdoms are not selected
  if ((which_taxonomy_rank+3) < 10) {
    
    unique_check<- unique_df_species[,(which_taxonomy_rank+3)] #save vector with +1 taxonomy rank than in which_taxonomy_rank (to_highlight and name still present in the dataframe)
    number_in_unique_check<- length(unique(unique_check)) #number of unique taxonomy ranks in taxonomy rank higher than which_taxonomy_rank
    
    #check if more selected elements in lowest taxonomy rank selection have the same higher taxonomy rank
    if (number_in_unique_check < 2) {
      
      lowest_selected_rank<- taxonomy_rank[which_taxonomy_rank+1]
      lowest_rank_name<- unique(unique_check)
      
    } else {
      taxonomy_checker<- FALSE
    }
    
  } else {
    taxonomy_checker<- FALSE
  }
}

if (taxonomy_checker == TRUE) {
  
  #this resaving is present due to the splitting of the script on two parts which in Shiny requires saving to reactive values (in TeloBase script named as values$). In the static version not required.
  values_lowest_selected_rank<- lowest_selected_rank
  values_lowest_rank_name<- lowest_rank_name
  values_unique_df_species<- unique_df_species
  
}

###

metacoder_sequence_to_plot<- "TTTAGGG" #sequence to highlight (in case of TeloBase it is based on the input) (part of the code is present due to the splitting of the script)

#check if some sequence to highlight was selected and change the colour for plotting
if(is.null(metacoder_sequence_to_plot)) {
  selected_metacoder_colour<- "#E6E6E3"
} else {
  selected_metacoder_colour<- "green" #colour for the highlighting (in case of TeloBase it is based on the input)
}

## This part of the code is present due to the splitting of the script on two parts -> one is present in the code related to the clicking of the button (above ###) and the second when rendering the plot itself 

  #using taxonomy filters to read the selected data for filtration of the unique_df_species (taxonomy filters were resolved in generate phylogenetic tree button code)
  
  lowest_selected_rank<- values_lowest_selected_rank
  lowest_rank_name<- values_lowest_rank_name
  
  #loading the dataframe as resolved in the generate phylogenetic tree button code
  
  unique_df_species<- values_unique_df_species
  
##

#taxonomy_rank vector to find position of lowest_selected_rank
taxonomy_rank<- c("species", "genus", "family", "order", "class", "phylum", "kingdom")

#position for lowest_rank_to_plot in taxonomy_rank
in_number_lowest_rank_to_plot<- which(taxonomy_rank == lowest_selected_rank)-2

#correct position if genus or species was selected so the in_number_lowest_rank_to_plot will be always at least 1
if(in_number_lowest_rank_to_plot < 1) {
  in_number_lowest_rank_to_plot<- 1
}

#translate in_number_lowest_rank_to_plot to the name by taxonomy_rank vector
lowest_rank_to_plot<- taxonomy_rank[in_number_lowest_rank_to_plot]

#filter unique_df_species to have only selected part based on the lowest_selected_rank
dataframe_lowest_rank_name<- unique_df_species[unique_df_species[[match(lowest_selected_rank, colnames(unique_df_species))]] == lowest_rank_name,]

#recalculate to_highlight based on lowest_rank_to_plot and further filter data to by the same vector
#which in lowest_rank_to_plot have the unique() to_highlight = 1
rank_to_highlight<- unique(dataframe_lowest_rank_name[[match(lowest_rank_to_plot, colnames(dataframe_lowest_rank_name))]][dataframe_lowest_rank_name$to_highlight == 1 & !is.na(dataframe_lowest_rank_name$to_highlight)])

dataframe_lowest_rank_name$to_highlight[dataframe_lowest_rank_name[[match(lowest_rank_to_plot, colnames(dataframe_lowest_rank_name))]] %in% rank_to_highlight]<- 1

unique_df_lowest_rank<- dataframe_lowest_rank_name[!duplicated(dataframe_lowest_rank_name[[match(lowest_rank_to_plot, colnames(dataframe_lowest_rank_name))]]),] #lower the amount of information for further plotting -> should help with the speed of the plotting

#filter NA sequences
unique_df_lowest_rank<- unique_df_lowest_rank[!is.na(unique_df_lowest_rank$to_highlight),]


#vector with taxonomy ranks to filter the point to which the plot will be done
taxonomy_rank<- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#filter taxonomy rank based on lowest_rank_to_plot
taxonomy_rank_filtered<- taxonomy_rank[1:match(lowest_rank_to_plot, taxonomy_rank)]

parsed_datafile<- parse_tax_data(unique_df_lowest_rank, class_cols = taxonomy_rank_filtered)

#sum the branches that needs to be highlighted and change the higher values to logical 0 and 1
logical_branches<- unlist(obs_apply(parsed_datafile, "tax_data", sum, value = "to_highlight"))
logical_branches[logical_branches > 0] <- 1

#connect to parsed_datafile new data set called summary_data that will be used for the highlight
parsed_datafile <- mutate_obs(parsed_datafile, data =  "summary_data",
                              taxon_id = taxon_ids,
                              branches_to_highlight = logical_branches
)

#heat_tree() plots randomly so to maintain a certain order, it is necessary to set a seed by set.seed() -> in several instances, the 123 does not have the best orientation of the plot as certain angles make the small graph
set.seed(123)

if(nrow(parsed_datafile$taxonomy_table()) == 1 & ncol(parsed_datafile$taxonomy_table()) == 3) {
  set.seed(8) 
}

if(nrow(parsed_datafile$taxonomy_table()) == 1 & ncol(parsed_datafile$taxonomy_table()) == 4) {
  set.seed(73) 
}

if(nrow(parsed_datafile$taxonomy_table()) == 1 & ncol(parsed_datafile$taxonomy_table()) == 5) {
  set.seed(64) 
}

if(nrow(parsed_datafile$taxonomy_table()) == 1 & ncol(parsed_datafile$taxonomy_table()) == 6) {
  set.seed(8) 
}

if(nrow(parsed_datafile$taxonomy_table()) == 1 & ncol(parsed_datafile$taxonomy_table()) == 7) {
  set.seed(88) 
}

#inputs in the TeloBase as (input$)
input_selected_metacoder_node_size<- 0.01
input_selected_metacoder_label_size<- 0.01
input_selected_metacoder_edge_size<- 0.001

heattreeplot<- heat_tree(parsed_datafile,
                         node_color = branches_to_highlight,
                         node_label = taxon_names,
                         initial_layout = "reingold-tilford", layout = "davidson-harel",
                         node_size_range = c(as.numeric(input_selected_metacoder_node_size), as.numeric(input_selected_metacoder_node_size)),
                         node_label_size_range = c(as.numeric(input_selected_metacoder_label_size), as.numeric(input_selected_metacoder_label_size)),
                         edge_size_range = c(as.numeric(input_selected_metacoder_edge_size), as.numeric(input_selected_metacoder_edge_size)),
                         node_color_range = c("#E6E6E3",selected_metacoder_colour),
                         node_color_interval = 0:1,
                         make_node_legend = FALSE
)

#create random name for temporarily stored .png file of the heattreeplot
random_number_name<- paste0(stri_rand_strings(1, 8), ".png")

#bind the name with the tempdir() path
filepath_temp_heattreeplot<- file.path(tempdir(),random_number_name)

#generate .png file
png(file = filepath_temp_heattreeplot, height = 3000, width = 3000)
print(heattreeplot)
dev.off()

#base64Encode the .png file for plotly output
txt <- RCurl::base64Encode(readBin(filepath_temp_heattreeplot, "raw", file.info(filepath_temp_heattreeplot)[1, "size"]), "txt")

#remove .png file from the tempdir()
file.remove(filepath_temp_heattreeplot)

#pixel size of the image (corresponds to generated png) and scale factor
image_width = 3000 
image_height = 3000 
scale_factor = 0.5 


#create empty plot_ly trace based on the image size parameters and scale_factor
metacoder_plot <- plot_ly() %>% 
  add_trace(x= c(0, image_width * scale_factor), 
            y= c(0, image_height * scale_factor), 
            type = 'scatter',  
            mode = 'markers', 
            alpha = 0)

#xy axis configuration
metacoder_plot<- metacoder_plot%>% layout(
  xaxis = list( 
    title = "", 
    zeroline = FALSE, 
    showline = FALSE, 
    showticklabels = FALSE, 
    showgrid = FALSE, 
    range = c(0, image_width * scale_factor) 
  ), 
  yaxis = list( 
    title = "", 
    zeroline = FALSE, 
    showline = FALSE, 
    showticklabels = FALSE, 
    showgrid = FALSE, 
    range = c(0, image_height * scale_factor), 
    scaleanchor="x" 
  )
) 

#add base encoded image to the plot
metacoder_plot<- metacoder_plot%>% layout( 
  images = list(  
    list(  
      source = paste('data:image/png;base64', txt, sep=','),  
      x=0, 
      sizex=image_width * scale_factor, 
      y=image_height * scale_factor, 
      sizey=image_height * scale_factor, 
      xref="x", 
      yref="y", 
      opacity=1.0, 
      layer="below", 
      sizing="stretch" 
    )  
  )) 

#configure the colour and plotly buttons
metacoder_plot<- config(metacoder_plot, displaylogo = FALSE, modeBarButtonsToRemove = c("toImage", "select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian")) %>% 
  layout(margin = list(left = 0, right = 0, bottom = 0, top = 0)) %>%
  layout(plot_bgcolor='#ffff',  
         xaxis = list(  
           zerolinecolor = '#ffff',  
           zerolinewidth = 2,  
           gridcolor = 'ffff'),  
         yaxis = list(  
           zerolinecolor = '#ffff',  
           zerolinewidth = 2,  
           gridcolor = 'ffff')  
  )

#show plot
metacoder_plot
