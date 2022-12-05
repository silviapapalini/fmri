library(tidyverse)

TR <- 2

subjects <- c("sub-FG01","sub-FG02","sub-FG03","sub-FG04","sub-FG05","sub-FG06","sub-FG07","sub-FG09","sub-FG10","sub-FG11","sub-FG12","sub-FG14","sub-FG16","sub-FG17","sub-FG20","sub-FG22","sub-FG23","sub-FG25","sub-FG28R","sub-FG29R","sub-FG30R","sub-FG31R","sub-FG32R","sub-FG33R","sub-RG26","sub-RG27","sub-RG30","sub-RG31","sub-RG34","sub-RG37","sub-RG38","sub-RG41","sub-RG44","sub-RG47","sub-RG48","sub-RG49","sub-RG50","sub-RG51","sub-RG55R","sub-RG56R","sub-RG58R","sub-RG59R","sub-RG61R","sub-RG62R","sub-RG63R","sub-RG64R")

read.subj <- function(subj) {
    events_file <- paste0(subj, "/func/", subj, "_task-Extinction_events.tsv")
    events <- read.delim(events_file, sep="\t") %>%
        dplyr::filter(trial_type != "relief_rating") %>%
        tidyr::separate(trial_type, into=c("timingcs", "cs"), sep="_") %>%
        dplyr::mutate(
            timingcs=factor(timingcs, labels=c("onset", "offset"), levels=c("onset", "relief")),
            cs=factor(cs, levels=c("csm", "csav", "csunav")))

    bold_file <- paste0("derivatives/afni/",subj,"/Extinction/", subj, "_task-Extinction_scaled_BOLD_from_NAcc.1D")
    bold <- unname(unlist(read.delim(bold_file, sep=" ", comment.char="#", header=FALSE)[1,-1]))

    fetch_window <- function(t) {
        start <- floor((t - 4) / TR)
        stop <- start + 20/TR
        x <- (bold[start:stop])
        names(x) <- seq(-4, by=TR, length.out=11)
        list(x)
    }
    X <- events %>%
        dplyr::rowwise() %>%
        dplyr::mutate(bold=fetch_window(onset)) %>%
        tidyr::unnest_longer(bold, indices_to="time")
    X
}
X <- lapply(subjects, read.subj)
names(X) <- subjects
X <- dplyr::bind_rows(X, .id="SUBJ_ID") %>%
    tidyr::extract(SUBJ_ID, "sub-([FR]G)([0-9]+R?)", into=c("group", "id"), remove=FALSE) %>%
    dplyr::mutate(early_late=ifelse(trial_number_cs <= 3, "early", ifelse(trial_number_cs <= 6, "middle", "late")))
write.table(X, file="onset_offset_bold.csv", sep=",", row.names=F, quote=F)

se <- function(x, na.rm=FALSE) sd(x, na.rm=na.rm) / sqrt(length(x))
S <- dplyr::group_by(X, group, timingcs, cs, time, early_late) %>% dplyr::summarize(mu=mean(bold, na.rm=TRUE), se=se(bold, na.rm=TRUE))
ggplot(S, aes(x=time, y=mu, color=cs, group=cs)) + geom_line() + geom_pointrange(aes(ymin=mu-se, ymax=mu+se), size=0.04) + facet_wrap(group+early_late ~ timingcs, scales="free_y")
