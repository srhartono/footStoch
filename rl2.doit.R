#nickPos = 1158 #1561
seed = 420
set.seed(seed)
nickPoss = c(0,679,703,1158,1561)
nickLengths = c(1,10,100,200)
for (nickPos in nickPoss) {
  for (nickLength in nickLengths) {
    #    nickPos = 703
    #    nickLength = 50
    minPeakLength = 50
    simN = 10000
    seqtypeswant.all = seqtypes
    seqtypeswant = c("pFC9")
    #seqtypeswant = c('FUS')
    #seqwant = myseqs[seqtypes == seqtypeswant]
    
    #all
    sigmas = c(-0.00,-0.01,-0.02,-0.03,-0.04,-0.05,-0.06,-0.07,-0.08,-0.09,-0.1,-0.14,seq(-0.00,-0.08,-0.002))
    #test
    # sigmas = c(-0.00,-0.01,-0.02,-0.03,-0.04,-0.05,-0.06,-0.07,-0.08,-0.09,-0.1,-0.14)
    # #679
    # sigmas = c(-0.00,-0.01,-0.02,-0.03,seq(-0.04,-0.05,by=-0.001),-0.06,-0.07,-0.08,-0.09,-0.1,-0.14)
    # #1158
    # sigmas = c(-0.00,-0.01,-0.02,seq(-0.03,-0.04,by=-0.001),-0.05,-0.06,-0.07,-0.08,-0.09,-0.1,-0.14)
    # #1561
    # sigmas = c(-0.00,-0.01,-0.02,-0.03,seq(-0.04,-0.05,by=-0.001),-0.06,-0.07,-0.08,-0.09,-0.1,-0.14)
    
    G.as = c(10)#0,4,8,10,15,20,25,50)#c(0,4,8,10,15,20,25,30)
    snames = gsub('^[ \\-]?0.(.+)$','\\1',ac(format(sigmas,digits=2))); snames
    
    #G.as = c(-10,0,5,10,15,20)#0,5,10,15,20,25,50,100)
    
    #grDevices::cairo_pdf(pasta(plasmid,'.pdf'),width=10,height=10)
    for (plasmid in seqtypeswant) {
      pdfname = pasta(plasmid)
      #if (file.exists(pdfname)) {
      #  print(paste(pdfname,"Nexted!"))
      #  next
      #}
      # next
      #      if (nickPos > 0) {
      pdfname = pasta(plasmid,'.nick',nickPos,'.nickLength',nickLength,'.',minPeakLength)
      #      }
      #pdf(pasta(pdfname,'.pdf'),width=10,height=50)
      pdf(pasta(pdfname,'_graphonly.pdf'),width=10,height=10)
      
      print(paste("Doing plasmid",plasmid))
      seq1 = myseqs[seqtypes == seqtypeswant]
      seq.len = nchar(seq1)
      strand = "+"
      start_from = NA
      
      bed = BEDS[BEDS$gene == plasmid,]
      bed = bed[grep('(Barcode|_diff|Amplicon)',bed$feature,invert=T),]
      if (grepl('T7_',plasmid)) {
        start_from = 'T7_Promoter'
        end_to     = 'RV_Primer'
        strand = '+'
      } else if (grepl('(pFC53|T3)',plasmid)) {
        start_from = 'T3_Promoter'
        end_to     = 'FW_Primer'
        strand = '-'
      } else {
        start_from = 'FW_Primer'
        end_to     = 'RV_Primer'
        strand = bed$strand[1]
        if (!defined(strand)) {
          strand = '+'
        }
        if (strand == '-') {
          end_to = 'FW_Primer'
          start_from = 'RV_Primer'
        }
      }
      
      if (strand == '-') {
        bed.end = max(G$n) - bed$beg
        bed.beg = max(G$n) - bed$end
        bed$beg = bed.beg
        bed$end = bed.end
      }      # G = G.a + G.B(i=n:n+m)) + G.S(m,sigma)
      
      ind = 0
      indtotal = length(G.as) * length(sigmas)
      # 1. Load sum(B)
      G.B = get_dgB(plasmid,myseqs,seqtypes)
      for (j in 1:length(G.as)) {
        for (s in 1:length(sigmas)) {
          ind = ind + 1
          G.a = G.as[j]
          sigma = sigmas[s]
          sname = snames[s]
          
          title0 = pasta('G of R-loop formation\nplasmid = ',plasmid)
          title0 = pasta(title0,' ; a = ',G.a)
          title0 = pasta(title0, ' ; sigma = ',sigma)
          if (nickPos > 0) {
            title0 = pasta(title0, ' ; nickPos = ',nickPos,' ; effectlength = ',nickLength)
          }
          title0 
          
          print(paste("Doing",ind,'/',indtotal,title0))
          name      = pasta('B.',sname,'.',G.a)
          name.top  = pasta(name,'.top')
          name.bot  = pasta(name,'.bot')
          name.prob = pasta(name,'.prob')
          name.plog = pasta(name,'.plog')
          
          saveDir = 'rl2/resources/dgsaves/'    
          saveFilename = pasta(saveDir,paste(plasmid,G.a,sname,nickPos,nickLength,minPeakLength,simN,sep='.'))
          saveFile.G.a_sigma = pasta(saveFilename,'.G.a_sigma.RDS')
          if (file.exists(saveFile.G.a_sigma)) {
            G.a_sigma = readRDS(saveFile.G.a_sigma)
          } else {
            G.a_sigma = calc.G_a_sigma(m.max = par$len.max, a = G.a, sigma = sigma, cons = cons, print = F);dim(G.a_sigma)[1]
            G.a_sigma = G.a_sigma[order(G.a_sigma$m),]
            saveRDS(G.a_sigma,saveFile.G.a_sigma)
          }    
          # Plot G.a_sigma: G vs rloop length at sigma, don't raelly need this but why not
          ylimcurr = -1 * 1300 #this is just to make it easy to flip plot
          if (ylimcurr < 0) {ylimmincurr = ylimcurr; ylimmaxcurr = 200} else {ylimmincurr = -200; ylimmaxcurr = ylimcurr}
          p0 = ggplot(G.a_sigma,aes(m, -1 * dg.sigma)) +
            geom_line(aes(m, -1 * dg.sigma)) +
            
            geom_line(aes(m, -1 * (dg.sigma+m*0.8)),color='red4') +
            annotate(geom='segment',x=0,xend=180,y=-1 * 1200,yend= -1 * 1200,color='red4') +
            annotate(geom='text',x=200,y=-1 * 1200,hjust=0,label='T Homopolymer',color='red4',size=3) +
            
            geom_line(aes(m, -1 * (dg.sigma+m*-0.1)),color='purple3') +
            annotate(geom='segment',x=0,xend=180,y=-1 * 1100,yend=-1 * 1100,color='purple3') +
            annotate(geom='text',x=200,y=-1 * 1100,hjust=0,label='Homogenized Random',color='purple3',size=3) +
            
            geom_line(aes(m, -1 * (dg.sigma+m*-0.2)),color='orange') +
            annotate(geom='segment',x=0,xend=180,y=-1 * 1000,yend=-1 * 1000,color='orange') +
            annotate(geom='text',x=200,y=-1 * 1000,hjust=0,label='Homogenized Favorable',color='orange',size=3) +
            
            geom_line(aes(m, -1 * (dg.sigma+m*-0.36)),color='green4') +
            annotate(geom='segment',x=0,xend=180,y=-1 * 900,yend=-1 * 900,color='green4') +
            annotate(geom='text',x=200,y=-1 * 900,hjust=0,label='G Homopolymer',color='green4',size=3) +
            
            coord_cartesian(xlim=c(0,1000),ylim=c(ylimmincurr,ylimmaxcurr)) +
            ggtitle(pasta('G of R-loop - sum(G.B), by R-loop length\na = ',G.a,' ; sigma = ',sigma)) +
            xlab('R-loop length (bp)') +
            ylab('-G (kcal/mol)') + theme_bw() + theme(panel.grid=element_blank()) +
            annotate(geom='segment',lty=2,x = 0,xend=1000,y= -1 *(G.a_sigma[G.a_sigma$m == 0,]$dg.sigma),yend= -1 * (G.a_sigma[G.a_sigma$m == 0,]$dg.sigma))

          if (ylimcurr < 0) {
            p0 = p0 + 
              annotate(geom='text',x = 1000,y=-1 * (G.a_sigma[G.a_sigma$m == 0,]$dg.sigma),vjust=0,hjust=1,label='above = forms R-loop',size=3)
          } else {
            p0 = p0 + 
              annotate(geom='text',x = 1000,y=-1 * (G.a_sigma[G.a_sigma$m == 0,]$dg.sigma),vjust=1,hjust=1,label='below = forms R-loop',size=3)
          }
          
          G = data.frame()
          my.xmin = 1
          my.xmax = 3000
          saveFile.G = pasta(saveFilename,'.G.RDS')
          if (file.exists(saveFile.G)) {
            G = readRDS(saveFile.G)
            my.xmax = max(G$n)-200
            
          } else {
            G   = subset(G.B,select=-B)
            rownames(G) = seq(1,dim(G)[1])
            G = G[order(G$n,G$m),]
            my.xmax = max(G$n)-200
            
            # - G.a = G.a[j]
            # - B(n,m) = G.B
            # - G.S = calc.G.a_sigma(sigma,G.a)
            G[,name]  = G.B$B + G.a_sigma$dg.sigma;dim(G)[1]
            
            #sanity check
            h(G)
            dim(G.B)[1] / dim(G.a_sigma)[1]
            dim(G)[1] / dim(G.a_sigma)[1]
            
            if (nickPos > 0) {
              G[G$m != 0 & G$n >= nickPos & G$n <= nickPos + nickLength,name] = G[G$m != 0 & G$n >= nickPos & G$n <= nickPos + nickLength,name] - G.a
            }
            # pos (n) = 1 to plasmid sequence length-500
            # len (m) = 0 to 2000
            # p(n.i,m.j) = e^(-G[pos.n.i,len.m.j]/RT) / sum(e^(-G[pos.n(1-max.i),len.m(1-max.j)]/RT))
            if (defined(start_from)) {
              if (defined(bed[bed$feature == start_from,])) {
                nstart = bed[bed$feature == start_from,]$beg
                G = G[G$n > nstart,]
              }
            }
            if (defined(end_to)) {
              if (defined(bed[bed$feature == end_to,])) {
                nend = bed[bed$feature == end_to,]$end - 100
                G = G[G$n < nend,]
              }
            }      
            G[,name.top]  = exp(-1 * G[,name]/(cons$R * cons$T)) # G of rloop at n=i, m=j
            G[,name.bot]  = sum(G[,name.top])                    # Sum of all G of rloop at n=(1-max.n), m=(1-max.m)
            G[,name.prob] = G[,name.top]/G[,name.bot]            # Probability of rloop at n=i, m=j
            h(G)
            # make it so probability of 0 becomes 4.9e-324
            #G[G[,name.prob] == 0,name.prob] = 4.940656e-324 #4.940656e-324lowest non zero R can do, which is smaller than .Machine$double.xmin (2.225074e-308)
            
            # turn probability into log for easier viewing
            # For probability=0, we can't do log(0) so use log (lowest non zero number R can do)
            # -> 4.940656e-324 is lowest non zero R can do, which is smaller than .Machine$double.xmin (2.225074e-308)
            # not using .Machine$double.xmin because lowest non zero R can do is smaller than that
            # (sometimes real probability is smaller than 2.2e-308 and using 2.2e-308 for p=0 will cause weird bumps)
            G[                  , name.plog] = 0
            G[G[,name.prob] != 0, name.plog] = log(G[G[,name.prob] != 0,name.prob],base=10)
            G[G[,name.prob] == 0, name.plog] = log(4.940656e-324, base = 10) 
            saveRDS(G,saveFile.G)
          }
          
          print("Graphing!")
          mwants = c(0,10,50,100,150,200,400,800,1000)
          G.mwant = G[G$m %in% mwants,]
          
          title0.small = gsub(' ','_',title0)
          title0.pdf = pasta(title0.small,'.pdf')
          
          # Graph energy
          G.mwant.y = -1 * G.mwant[,name]
          G.mwant.y.max = max(G.mwant.y)#max(abs(G.mwant.y)*0.1) 
          G.mwant.y.min = -400 #min(G.mwant[G.mwant$m == 0,name])-100
          
          p.G = ggplot(G.mwant,aes(x=n,y=1)) +
            theme_bw() + theme(panel.grid=element_blank(),legend.position='bottom',legend.box = 'vertical') +
            guides(color = guide_legend(title='R-loop Length',order = 1),fill='none') +
            ggtitle(title0) + ylab('-G (kcal/mol)') + xlab("Pos in plasmid (bp)") +
            geom_line(aes(y=G.mwant.y,color=af(m))) + 
            coord_cartesian(xlim=c(1,my.xmax),ylim=c(G.mwant.y.min,200)) +
            geom_rect(data=bed,aes(x=beg,y=0,xmin=beg,xmax=end,ymin=G.mwant.y.min,ymax = 200,fill=feature),alpha=0.15,color=NA) +
            geom_text(data=bed,aes(x=(beg+end)/2,y=200,label=feature),angle=45,size=2,hjust=1,color='black')
          
          
          # graph probability
          prob.mwant.y = G.mwant[,name.prob]
          prob.mwant.y.max = max(prob.mwant.y)# + max(abs(prob.mwant.y)*0.2) 

          p.prob = ggplot(G.mwant,aes(x=n,y=1)) +
            theme_bw() + theme(panel.grid=element_blank(),legend.position='bottom',legend.box = 'vertical') +
            guides(color = guide_legend(title='R-loop Length',order = 1),fill='none') +
            ggtitle(pasta('Prob of R-loop formation\na = ',G.a,' ; sigma = ',sigma)) +
            ylab('Probability') + xlab("Pos in plasmid (bp)") +
            geom_line(aes(y=G.mwant[,name.prob],color=af(m))) +
            geom_line(aes(y=prob.mwant.y,color=af(m))) + 
            coord_cartesian(xlim=c(1,my.xmax),ylim=c(min(prob.mwant.y),max(prob.mwant.y.max))) +
            geom_rect(data=bed,aes(x=beg,y=0,xmin=beg,xmax=end,ymin=min(prob.mwant.y),ymax = (prob.mwant.y.max),fill=feature),alpha=0.15,color=NA) +
            geom_text(data=bed,aes(x=(beg+end)/2,y=max(prob.mwant.y),label=feature),angle=45,size=2,hjust=1,color='black')
          
          # Simulation
          saveFile.df.sim = pasta(saveFilename,'.df.sim.RDS')
          if (file.exists(saveFile.df.sim)) {
            df.sim = readRDS(saveFile.df.sim)
          } else {
            df.sim = G[sample(x=seq(1,dim(G)[1]),simN,prob = G[,name.prob],replace=T),]
            saveRDS(df.sim,saveFile.df.sim)
          }    
          
          rloop.len = paste('Rloop len at 0%, 25%, 50%, 75%, and 100%: ',paste(quantile(df.sim$end-df.sim$beg,probs=c(0,0.25,0.5,0.75,1)),collapse=','))
          perc1 = ai(dim(df.sim[df.sim$m > minPeakLength,])[1] / simN * 1000+0.5)/10
          #      title1 = paste('plasmid',plasmid,',σ,a=',G.a)
          if (dim(df.sim[df.sim$m > minPeakLength,])[1] > 0) {
            df.sim = df.sim[df.sim$m > minPeakLength,]
            df.sim$beg = df.sim$n
            df.sim$end = df.sim$beg + df.sim$m
            
            df.sim = df.sim[order(df.sim$beg,df.sim$end),]
            df.sim$y = seq(1,dim(df.sim)[1])
            yused3 = df.sim$y
            yused3.max = max(yused3)
            
            p1len = ggplot(df.sim,aes(end-beg)) +
              geom_boxplot(outlier.shape=NA) +
              theme_bw() + theme(panel.grid=element_blank()) +
              annotate(geom='text',x=median(df.sim$end-df.sim$beg),y=0,hjust=0,size=7,label=median(df.sim$end-df.sim$beg)) +
              ggtitle(rloop.len)
            
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
            myPalette2 = rep(myPalette(5000),1)
            #          scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))
            p1 = ggplot(df.sim,aes(beg,y)) + geom_segment(aes(x=beg,xend=end,y=y,yend=y,color=end-beg)) + 
              theme_bw() + 
              theme(panel.grid=element_blank(),axis.text.x=element_text(angle=90),legend.position='bottom',legend.box = 'vertical') +
              guides(color = guide_legend(title='R-loop Length',order = 1),fill='none') +
              scale_x_continuous(breaks = seq(0,my.xmax,200)) +
              coord_cartesian(xlim = c(1,my.xmax),ylim=c(0,yused3.max)) +
              geom_rect(data=bed,aes(x=beg,y=0,xmin=beg,xmax=end,ymin=min(yused3),ymax = yused3.max,fill=feature),alpha=0.15,color=NA) +
              geom_text(data=bed,aes(x=(beg+end)/2,y=yused3.max,label=feature),angle=45,size=2,hjust=1,color='black') +
              scale_color_gradientn(colors=myPalette2,limits=c(1,5000))
            p1 = p1 + ggtitle(pasta('10k Simulated Reads (',perc1,' %peak )\n',rloop.len,'\nplasmid = ',plasmid,' ; a = ',G.a,' ; sigma = ',sigma,' ; nickPos = ',nickPos,' ; effectlength = ',nickLength)) + 
              ylab('Read Number (sorted by Rloop starts)') + xlab("Pos in plasmid (bp)")

          } else {
            
            testfail = data.frame(beg=seq(1,my.xmax),end=seq(1,my.xmax))
            p1len = ggplot(testfail,aes(end-beg)) +
              geom_boxplot(outlier.shape=NA) +
              theme_bw() + theme(panel.grid=element_blank()) +
              ggtitle(paste('No rloop'))
            
            
            p1 = ggplot(testfail,aes(beg,end)) +
              theme_bw() + theme(panel.grid=element_blank(),legend.position='bottom',legend.box = 'vertical') +
              guides(color = guide_legend(title='R-loop Length',order = 1),fill='none') +
              coord_cartesian(xlim = c(1,my.xmax),ylim=c(0,my.xmax)) +
              geom_line() +
              annotate(geom='text',x=my.xmax/2,y=my.xmax/2,label='No Rloop very sad') +
              geom_rect(data=bed,aes(x=beg,y=0,xmin=beg,xmax=end,ymin=1,ymax =my.xmax,fill=feature),alpha=0.15,color=NA) +
              geom_text(data=bed,aes(x=(beg+end)/2,y=my.xmax,label=feature),angle=45,size=2,hjust=1,color='black') +
              ggtitle(pasta('10k Simulated Reads\nplasmid = ',plasmid,' ; a = ',G.a,' ; sigma = ',sigma)) + ylab('Read Number (sorted by Rloop starts)') + xlab("Pos in plasmid (bp)")
            # p2 = ggplot(data.frame(x=seq(1,10),y=seq(1,10)),aes(x,y)) + geom_text(x=5,y=5,label='No Rloop very sad')
          }
          # grid.arrange(p1,p1len,p.G,p.prob,p0,ncol=1)
          grid.arrange(p1)
        }
      }
      dev.off()
    }
  }
}