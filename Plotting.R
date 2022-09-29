

data_op <- as.data.frame(cbind(pred_j, validation$Age))

p2 <- ggplot(data_op,aes(x =  V2, y = pred_j))+geom_point() +
  xlim(0,80) + ylim(0,80)+ geom_abline(slope=1, intercept=0, linetype = 2) +
  annotate("text", x = 60, y = 20, label= "RMSE = 6.63")+          
  annotate("text", x = 60, y = 10, label= "MAD = 4.88")

p2 + xlab("Chronological Age (years)") + ylab("Predicted Age (Years)") + theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Whole Blood chrX + autosomal probes, Male, optimal model")

