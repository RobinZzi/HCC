# 加载必要的包
library(ggplot2)
library(dplyr)
rm(list=ls())
save.image("ga_methcomb.Rdata")
# 生成模拟数据（实际使用时替换为你的数据框）
test_ga <-test
# 计算象限百分比
test_ga <- test_ga %>%
  mutate(quadrant = case_when(
    RNAlog2FC > 0 & DNA_meth_dt > 0 ~ "Q1",
    RNAlog2FC > 0 & DNA_meth_dt < 0 ~ "Q2",
    RNAlog2FC < 0 & DNA_meth_dt < 0 ~ "Q3",
    RNAlog2FC < 0 & DNA_meth_dt > 0 ~ "Q4",
    TRUE ~ "Axis"
  ))

test_ga <- subset(test_ga,subset=quadrant!="Axis")


percentages <- test_ga %>%
  group_by(quadrant) %>%
  summarise(n = n()) %>%
  mutate(percent = round(100 * n / sum(n), 1))


xmax <- 1000
xmin <- -1000
ymax <- max(test_ga$DNA_meth_dt)
ymin <- min(test_ga$DNA_meth_dt)

# 定义象限标注位置
label_positions <- data.frame(
  quadrant = c("Q1", "Q2", "Q3", "Q4"),
  y = c(xmax*0.95, xmin*0.95, xmin*0.95, xmax*0.95),
  x = c(ymax*0.95, ymax*0.95, ymin*0.95, ymin*0.95)
) %>%
  left_join(percentages, by = "quadrant") %>%
  mutate(label = paste0(quadrant, ": ", percent, "%"))



ggplot(test_ga, aes(y = RNAlog2FC, x = DNA_meth_dt,color=quadrant)) +
  geom_point(alpha = 0.4, size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw()+xlim(-1,1)+scale_color_brewer(palette = "Set2")
