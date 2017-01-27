####################################################
# script by Jonathan Rolland and Daniele Silvestro #
####################################################

"This script containt the functions and algorithms implemented to estimate ancestral states under a
Brownian motion model of evolution. It uses phylogenies of extant taxa and trait values at the tips
as well as informative priors on the trait values at internal nodes derived from fossil data.

Full documentation of the script is still in development but it can be used using the example files
provided here https://github.com/jonathanrolland/niche_evolution. To run the script change the settings below and then execute the entire
script in R.

The script requires R and the R libraries 'phytools','diversitree','gtools'."

require(phytools)
require(diversitree)
library(gtools)

######################################################################################################################################################
#########################################		         RUN THE MODELS			  	##################################################################
######################################################################################################################################################

w_dir 			= getwd() # set working directory for data
wdir2			= getwd() # set working directory for output

#choose your group between 1: birds, 2: mammals, 3: amphibians, 4: squamata.
group<-1

#choose the variable you want to reconstruct
col_choose<-1
col_distrib_data<-c("alt_min","alt_max","latitude_min","latitude_max","min_temp","max_temp")[col_choose]

# load the empirical data
# to run this correctly, open R in the data directory.

if(group==1){
tree_file        = "birds.tre"
fossil_tbl_file  = "fossil_birds.txt"
distr_tbl_file   = "birds_data.txt"
fossil_file_tbl  = "fossil_tbl_birds_save.txt"
}

if(group==2){
tree_file        = "mammals.tre"
fossil_tbl_file  = "fossil.mammals.txt"
distr_tbl_file   = "mammals_data.txt"
fossil_file_tbl  = "fossil_tbl_mammals_save.txt"
}

if(group==3){
    tree_file        = "amphibians.tre"
    fossil_tbl_file  = "fossil_amphibia.txt"
    distr_tbl_file   = "amphibians_data.txt"
    fossil_file_tbl  = "fossil_tbl_amphibians_save.txt"
}

if(group==4){
    tree_file        = "squamates.tre"
    fossil_tbl_file  = "fossil_squamata.txt"
    distr_tbl_file   = "squamates_data.txt"
    fossil_file_tbl  = "fossil_tbl_squamata_save.txt"
}

#parameters to run the models

nGibbs_gen      = 10000000 # total number of mcmc generations 
Gibbs_sample    = 1000 # sampling frequency 
print_f         = 1000000 # print frequency in terminal

#sample			= 100 # record samples every xx generations

######################################################################################################################################################
#########################################		         FUNCTIONS 			    	##################################################################
######################################################################################################################################################

dynamicPlot     = F # functionnality with plot (do not use here)
useHP           = F # functionnality with Hyper prior (do not use here)

update_multiplier_proposal <- function(i,d){
    u = runif(1,0,1)
    l = 2*log(d)
    m = exp(l*(u-0.5))
    ii = i * m
    U=log(m)
    return(c(ii, U))
}

fN<-function(xa, xb, vpa, vpb, sigma2) {
    #the normal density. Sigma2 is the variance. Same as log(dnorm(xa-xb, 0, sqrt(sigma * sigma * (vpa + vpb))))
    return( dnorm((xa-xb),mean=0,sd= sqrt(sigma2 * (vpa + vpb)), log=T) )
    #return((-(xa - xb)^2 / (2 * sigma2 * (vpa + vpb)) - log(sqrt(2 * pi * sigma2 * (vpa + vpb)))));
}

fN2<-function(xa, xb, vpa, vpb, sigma2_a,sigma2_b, anc) {
    #same as fN but with a known ancestor instead of xa-xb, see Joe's book eqn 23.10 (some mistakes in the denominator though in his equation)
    return(dnorm(xa, mean=anc, sd= sqrt(vpa*sigma2_a),log=T) + dnorm(xb, mean=anc, sd= sqrt(vpb*sigma2_b),log=T)  )
    #return(-(((xa - anc)^2 / vpa) + ((xb - anc)^2 / vpb)) / (2 * sigma2) - log(sqrt(2 * pi * sigma2 * (vpa + vpb))));
}

newlnLike<-function(tree, vector_tip_root_nodes_values, sigma2,D) {

    vec_values = vector_tip_root_nodes_values
    anc_ind_v <- D[,1]
    a_ind_v   <- D[,2]
    b_ind_v   <- D[,3]
    vpa_v     <- D[,4]
    vpb_v     <- D[,5]
    
    anc_v <- vec_values[anc_ind_v]
    a_v   <- vec_values[a_ind_v]
    b_v   <- vec_values[b_ind_v]
    
    a_ind_v[a_ind_v>ntips] =a_ind_v[a_ind_v>ntips]-1
    b_ind_v[b_ind_v>ntips] =b_ind_v[b_ind_v>ntips]-1
    
    L1 = dnorm(a_v, mean=anc_v, sd= sqrt(vpa_v*sigma2[a_ind_v]),log=T)
    L2 = dnorm(b_v, mean=anc_v, sd= sqrt(vpb_v*sigma2[b_ind_v]),log=T)
    
    L = L1+L2
    return(L)
    
}

build_table <- function (tree,ntips,data){
    
    tree.edge.ordered<-tree$edge[order(tree$edge[,2]),]
    # the root is at ntips+1
    root<-ntips+1
    
    # the table with ancestor, the two descendants and the two branch length
    table.tree<-matrix(NA,ntips-1,5)
    colnames(table.tree)<-c("ancestor","descendant1","descendant2","branch.length1", "branch.length2")
    
    #contains the root and all the nodes
    table.tree[,1]<-c(ntips+1,tree.edge.ordered[(ntips+1):(ntips+ntips-2),2])
    
    for(i in 1:dim(table.tree)[1]){
        table.tree[i,c("descendant1","descendant2")] <- tree$edge[which(tree$edge[,1]==table.tree[i,1]),2]
        table.tree[i,c("branch.length1","branch.length2")] <- tree$edge.length[which(tree$edge[,1]==table.tree[i,1])]
    }
    return(table.tree)
}

get_calibration_tbl <- function (D,root_calibration){
    calibration_tbl =NULL
    for (i in 1:dim(D)[1]){
        calibration_tbl = rbind(calibration_tbl,c(0,1))
    }
    row.names(calibration_tbl)=D[,1]
    calibration_tbl[1,] = root_calibration
    return(calibration_tbl)
}

index_prior <- function(tree){
    bt<-branching.times(tree)
    names(bt)<-(length(tree$tip.label)+1):(length(tree$tip.label)*2-1)
    
    table<-matrix(0,length(bt),2)
    colnames(table)<-c("node","age")
    table[,1]<-names(bt)
    table[,2]<-bt
    
    table2<-matrix(0,length(bt),2)
    colnames(table2)<-c("prior","category")
    
    table2[,1]<-0
    table2<-(apply(table2,2,as.numeric))
    
    #time slice
    for (i in 1:dim(table2)[1]){
        node_age = as.numeric(table[i,2])
        # GEOLOGIC TIME SCALE 2013
        if (node_age< 2588){
            table2[i,2] = 1
        }else if  (node_age<23){
            table2[i,2] = 2
        }else if (node_age<66){
            table2[i,2] = 3
        }else if (node_age<145){
            table2[i,2] = 4
        }else if (node_age<201.3){
            table2[i,2] = 5
        }else if (node_age<252.17){
            table2[i,2] = 6
        }else{
            table2[i,2] = 7
        }
    }
    
    rownames(table2)<- table[,1]
    return (table2)
    
}

get_fossil_table <- function(tree,fossil_tbl_file, distr_tbl_file){
	calc_mean_sd_species <- function(distr_tbl_file){
		distrib_raw_data = read.csv2(distr_tbl_file,head=T,sep="\t")
		M = distrib_raw_data[,"latitude_max"]
		m = distrib_raw_data[,"latitude_min"]
		vector <- c(M-m)
		mean_sd = mean(vector)/(2*1.96) 		
		# meanrange = mean(M-m) = mu + 1.96*sigma/sqrt(n) - (mu - 1.96*sigma/sqrt(n))= 2*1.96sigma/sqrt(n)
		# because n=1; sigma = mean(latmin - latmax)/2*1.96
		
		return(mean_sd)
	}

	tiplab = tree$tip.label
	genera_list_tree = rep(0,length(tiplab))
	for (i in 1:length(tiplab)){
		genera_list_tree[i] = strsplit(tiplab[i],"_")[[1]][1]
	}	

	fossil_raw_data = read.table(fossil_tbl_file,header = T,sep="\t")
	genera_list_fossil= unique(fossil_raw_data[,"occurrence.genus_name"])
	genera_list_fossil = as.character(genera_list_fossil)
	mean_sd_species = calc_mean_sd_species(distr_tbl_file)

	fossil_tbl = NULL
	BT<-branching.times(tree)
	
	for ( i in 1:length(genera_list_fossil)){
		genus = genera_list_fossil[i]
		genus_ind =which(fossil_raw_data[,"occurrence.genus_name"]==genus )
		fossil_occs_age = fossil_raw_data[genus_ind,"ma_mid"]
		fossil_occs_lat = fossil_raw_data[genus_ind,"paleolatdec"]
	
		phylo_genus_tips = tiplab[which(genera_list_tree == genus)]
	
		if (length(phylo_genus_tips)>0){
			if (length(phylo_genus_tips)==1){
				phylo_genus_tips = c(phylo_genus_tips,phylo_genus_tips)
				list<-as.numeric(names(table(findMRCA(tree, tips= phylo_genus_tips, type="node"))))
				genus_mrca_node = tree$edge[order(tree$edge[,2]),][list,1]
			}else{
				list<-as.numeric(names(table(findMRCA(tree, tips= phylo_genus_tips, type="node"))))
				 genus_mrca_node = min(list[which(list>length(tree$tip.label))])
			 }
				
			genus_mrca_age   = BT[genus_mrca_node-length(tree$tip.label)]
			
			delta_t = abs(fossil_occs_age-genus_mrca_age)
			WMean = sum(fossil_occs_lat*delta_t)/sum(delta_t)
			E = matrix(c(genus_mrca_node,WMean, mean_sd_species),1,3)
			fossil_tbl = rbind(fossil_tbl, E)
		}
	
	}
	return(fossil_tbl)
}

################ RUN GIBBS FUNCTIONS
get_joint_mu_s2 <- function (mu_f,s2_f,mu_g,s2_g){
    s2_fg = (s2_f*s2_g)/(s2_f+s2_g)
    mu_fg = (mu_f/s2_f + mu_g/s2_g) * s2_fg
    return(c(mu_fg,s2_fg))
}

runGibbs <- function(sigma2, vector_tip_root_nodes_values,D,prior_tbl,get_expected=0) {
    vec_values = vector_tip_root_nodes_values
    # loop over the table from most recent to root
    for (i in dim(D)[1]:2){
        #for (rep in 1:20){
        #	i = sample(dim(D)[1])[1]
	      	anc_ind <- D[i,1];
            a_ind   <- D[i,2];  # index of descendants
            b_ind   <- D[i,3];  # index of descendants
            vpa     <- D[i,4];  # br length
            vpb     <- D[i,5];  # br length
            anc = vec_values[anc_ind]
            a = vec_values[a_ind]
            b = vec_values[b_ind]
            # calibration prior
            calibrated_prior_mu = prior_tbl[which(rownames(prior_tbl)==anc_ind),1]
            calibrated_prior_s2 = prior_tbl[which(rownames(prior_tbl)==anc_ind),2]
            
            if (a_ind>ntips){
                sig_a_ind=a_ind-1 # because no sig2 value for root
            }else{sig_a_ind=a_ind}
            
            if (b_ind>ntips){
                sig_b_ind=b_ind-1
            }else{sig_b_ind=b_ind}
            
            desc_lik = get_joint_mu_s2(a, (vpa)*sigma2[sig_a_ind], b, (vpb)*sigma2[sig_b_ind])
            prior_prm = prior_tbl[which(rownames(prior_tbl)==anc_ind),]
            
            # lik from (stem) ancestral state
            if (i>1){
                if (anc_ind %in% D[,2]){
                    stem_ind = D[which(D[,2]==anc_ind),1]
                    stem_brl = D[which(D[,2]==anc_ind),4]
                }
                else{
                    stem_ind = D[which(D[,3]==anc_ind),1]
                    stem_brl = D[which(D[,3]==anc_ind),5]
                }
                
                sig_stem_ind=stem_ind-1

                stem_val = vec_values[stem_ind]
                stem_lik = c(stem_val, (stem_brl)*sigma2[sig_stem_ind])
                
                lik_val = get_joint_mu_s2(stem_lik[1],stem_lik[2],desc_lik[1],desc_lik[2])
                
                post_prm= get_joint_mu_s2(lik_val[1],lik_val[2],prior_prm[1], prior_prm[2])
                if (get_expected==1){
                    vec_values[anc_ind] = post_prm[1]
                }else{
                    vec_values[anc_ind] = rnorm(1, mean=post_prm[1], sd=sqrt(post_prm[2]))
                }
            }
    }
    return(vec_values[D[,1]] )
}


################################## START MCMC ###################################
mcmc.gibbs4 <- function (tree, x, alter_ind, D, prior_tbl_ind, fossil_tbl, true_anc= 0, ngen = 100000, control = list(), gibbs_sampler=T, useVCV=F, sample=sample, logfile="log", dynamicPlot = F, useFixPart=F, useHP=T){
    
    time_frames = max(prior_tbl_ind[,2])
    if (useHP==F){hp_sd = rep(10000, time_frames)
    }else{hp_sd = rep(10, time_frames)}
    hp_mu = rep(0, time_frames)
    
    cat(c("it", "posterior","likelihood","prior", "sig2",paste("mu_hp", 1:time_frames ,sep="_"),paste("sig_hp", 1:time_frames ,sep="_"),"root",D[-1,1], "\n"),sep="\t", file=logfile, append=F) # , paste("sig2", 1:length(TE[,1]),sep="_"), D[-1,1],"\n")
    # starting value for root: mean of all x values
    a <- mean(x)
    # starting values for anc states: random within the central 33% of x values
    x_temp = sort(x)
    y <- runif(tree$Nnode - 1, x_temp[as.integer(0.33*length(x))],x_temp[as.integer(0.66*length(x))])
    # starting values sig2: st dev of x values
    sig2 <- sd(x)/3
    
    ind_sig2 <- rep(1, length(c(x, y))) #1:length(c(x, y)) #
    mean_sig <- rep(0, length(c(x, y)))
    mean_anc <- rep(0, length(c(a, y)))
    
    # lik function
    if (useVCV ==T){
        temp <- phyl.vcv(as.matrix(x), vcv(tree), 1)
        likelihood <- function(C, invC, detC, x, sig2, a, y) {
            z <- c(x, y) - a
            logLik <- -z %*% invC %*% z/(2 * sig2) - nrow(C) * log(2 * pi)/2 - nrow(C) * log(sig2)/2 - detC/2
            return(logLik)
        }
        C <- vcvPhylo(tree)
        if (any(tree$edge.length <= (10 * .Machine$double.eps)))
        stop("some branch lengths are 0 or nearly zero")
        invC <- solve(C)
        detC <- determinant(C, logarithm = TRUE)$modulus[1]
    }else{
        likelihood <- function(C, invC, detC, x, sig2, a, y) {
            vector_tip_root_nodes_values = c(x, a, y)
            return(sum(newlnLike(tree, vector_tip_root_nodes_values, sig2,D)))
        }
    }
    
    prior_tbl=prior_tbl_ind
    prior_tbl[,2] = hp_sd[prior_tbl_ind[,2]]
    if (!is.null(fossil_tbl) ){prior_tbl[rownames(fossil_tbl),] = fossil_tbl[,1:2]}
    print(head(prior_tbl))
	
	
    # priors
    log.prior <- function(sig2, a, y, prior_tbl) {
        prior_sig2 <- sum(dcauchy(sig2, loc=0,scale=1, log = TRUE) )
        prior_anc = sum(dnorm(c(a, y), mean = prior_tbl[,1], sd = prior_tbl[,2], log = T))
        return(prior_sig2+prior_anc)
    }
    
    x <- x[tree$tip.label]
    if (is.null(names(y))) {
        names(y) <- length(tree$tip) + 2:tree$Nnode
    }else {y[as.character(length(tree$tip) + 2:tree$Nnode)]}
    
    L  <- newlnLike(tree, c(x, a, y), sig2[ind_sig2],D)
    Pr <- log.prior(sig2, a, y,prior_tbl)
    
    
    IND_edge = c()
    for (ind_edge in tree$edge[,2]){
        if (ind_edge>ntips){
            ind_edge_t=ind_edge-1
        }else{ind_edge_t=ind_edge}
        IND_edge = append(IND_edge,ind_edge_t)
    }
    
    
    # START MCMC
    for (i in 1:ngen) {
        y.prime     = y
        L.prime     = L
        Pr.prime    = Pr
        a.prime     = a
        sig2.prime  = sig2
        gibbs=0
        hastings=0
        
        if (i%%print_f == 0 || i==1) {
            cat(c("\n",i,round(c(sum(L),Pr,sig2),2))) #,"\n",ind_sig2 ))
        }
        
        j <- (i - 1)%%(tree$Nnode + 1)
        rr= runif(1,0,1)
        
        update_sig_freq=0.85
        if (rr<update_sig_freq) {
            if (runif(1,0,1)>0.5){
                if (runif(1,0,1)<0.5){
                    sig2_update <-  update_multiplier_proposal(sig2.prime,1.1)
                    sig2.prime = sig2_update[1]
                    hastings = sig2_update[2]
                }else{
					sig2.prime = abs(sig2 + rnorm(n = 1, sd = sample(c(0.5,1,2,5),1)  ))
                    
                }
            }
            else{
                a.prime <- a + rnorm(n = 1, sd = (0.5*sd(x)/3)) #sqrt(con$prop[j + 1]))
            }
        }
        else{
            # ANC STATES
            if (gibbs_sampler==F){
                k <- j - 1
                y.prime <- y
                y.prime[k] <- y[k] + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))
            }
            else {
                vector_tip_root_nodes_values = c(x, a.prime, y.prime)

                    vector_tip_root_nodes_values = c(x, a.prime, y.prime)
                    
                    y_temp = runGibbs(sig2.prime[ind_sig2], vector_tip_root_nodes_values,D,prior_tbl,get_expected=0)
					print(head(prior_tbl))
					print(head(x))
                    a.prime= y_temp[1]
                    y.prime= y_temp[-1]
                gibbs=1
            }
            
        }
        
	       # calc post
           L.prime <- newlnLike(tree, c(x, a.prime, y.prime), sig2.prime[ind_sig2],D)
           Pr.prime <- log.prior(sig2.prime, a.prime, y.prime, prior_tbl)
           
           if ( (sum(Pr.prime) + sum(L.prime) - sum(Pr) - sum(L) + hastings) >= log(runif(1,0,1)) || gibbs==1){
               y =    y.prime
               L =    L.prime
               Pr =   Pr.prime
               a =    a.prime
               sig2 = sig2.prime
           }
           
           if (i%%sample == 0) {
               rates_temp=sig2[ind_sig2]
               # rel_err =  (c(a, y) - true_anc)
               #MAE = mean(abs(c(a, y)-true_anc))
               MAPE = sum(abs((rates_temp[IND_edge]-true_anc)))
               cat(c(i,sum(L)+sum(Pr), sum(L),sum(Pr), mean(sig2[ind_sig2]), hp_mu, hp_sd, a, y, "\n"),sep="\t", file=logfile, append=T)
           }
    }
}

##############
#set directory
setwd(w_dir)

#load the trees
if(group==2) tree = read.nexus(tree_file)
if(group==1|group==3|group==4) tree = read.tree(tree_file)

#load the data of the distribution for each species
raw_data= read.csv2(distr_tbl_file,head=T,sep="\t")

#Renames the species to uniformize tables
Name_sp<-rep(NA,dim(raw_data)[1])
for( i in 1: dim(raw_data)[1]){
	z<-strsplit(as.character(raw_data$Name_sp[i])," ")
	Name_sp[i]<-as.character(paste(z[[1]][1],"_",z[[1]][2],sep=""))
}

#delete all the tips of the tree without species name (check)
tree<-drop.tip(tree,tree$tip.label[which(tree$tip.label%in%Name_sp==F)])

#check on the labels of the nodes
tree$node.label<-c((length(tree$tip.label)+1): (length(tree$tip.label)*2-1))

#count the number of tips
ntips = length(tree$tip.label)

raw_data[,1] = Name_sp
data = raw_data[,col_distrib_data] 

#for altitude we divided by 1000 for facilitating the estimation of sigma. 
if(col_choose<3 ){ data = raw_data[,col_distrib_data]/1000}
names(data) = raw_data[,"Name_sp"]

D <- build_table(tree, ntips,data)
prior_tbl_ind = index_prior(tree)
prior_tbl_ind[,1]= mean(data)

#in the cases where you reconstruct traits that needs fossil data (e.g. latitude)

### Uncomment these lines if it is the first time that you run the code (this will open the raw files containing the fossil information). 
### This will calculate and write a fossil table : fossil_tbl summarizing the priors that will then be used in the reconstruction.
### As this code is computer intensive, the best is to compute the table only once (the output files are provided here for the four groups: "fossil_tbl_birds_save.txt",...).
#fossil_tbl  = get_fossil_table(tree,fossil_tbl_file, distr_tbl_file)
#print(fossil_tbl)
#colnames(fossil_tbl) = c("node","mean","sd")
#write.table(fossil_tbl, file=fossil_file_tbl,sep="\t")

# uncomment this if you have already computed the fossil_tbl
fossil_tbl = NULL
if(col_choose==3 |col_choose==4){
	fossil_tbl_temp = read.table(fossil_file_tbl,h=T,sep="\t")
#use the information contained in fossil_tbl_temp to create a table in the right format to use mcmc function
for (i in unique(fossil_tbl_temp$node)){
	row = fossil_tbl_temp[fossil_tbl_temp$node==i,]
	row=apply(row,2,FUN=mean)
	row[1] = as.integer(row[1])
	fossil_tbl = rbind(row, fossil_tbl)}
rownames(fossil_tbl) = fossil_tbl[,1]
fossil_tbl = fossil_tbl[,-1]
}

#change directory for the output
setwd(wdir2)

#create a new logfile that will contain the results of the mcmc
logfile3 = sprintf("group%smcmc_test.log", paste(group,col_distrib_data,sep=""))

#run the mcmc
mcmc.gibbs4(tree, data, alter_ind, D, prior_tbl_ind, fossil_tbl, true_anc=0, ngen= nGibbs_gen, gibbs_sampler=T, logfile=logfile3, sample=Gibbs_sample)


