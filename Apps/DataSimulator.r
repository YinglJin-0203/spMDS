# ============================================================
# Longitudinal Multivariate Dataset Generator
# User-specified number of variables x time points x subjects
# Pairwise correlations drift sinusoidally over time.
# ============================================================
#
# Required packages:
#   install.packages(c("shiny", "ggplot2", "tidyr", "dplyr", "DT", "Matrix"))

library(shiny)
library(ggplot2)
library(tidyr)
library(dplyr)
library(DT)
library(Matrix)

VAR_COLORS_PALETTE <- c(
  "#6366f1", "#ec4899", "#10b981", "#f59e0b", "#06b6d4",
  "#a855f7", "#ef4444", "#84cc16", "#3b82f6", "#f97316",
  "#14b8a6", "#eab308", "#8b5cf6", "#f43f5e", "#22c55e"
)

var_names_for <- function(n_vars) paste0("V", 1:n_vars)
var_color_for <- function(i, n_vars) VAR_COLORS_PALETTE[((i - 1) %% length(VAR_COLORS_PALETTE)) + 1]

pairs_for <- function(n_vars) combn(n_vars, 2, simplify = FALSE)
pair_label <- function(p, var_names) paste0(var_names[p[1]], "-", var_names[p[2]])

# ---------- Nearest valid correlation matrix ----------
# If sinusoidally-drifting pairwise correlations combine into a
# matrix that isn't positive semi-definite, project it to the
# nearest valid correlation matrix (clip negative eigenvalues,
# then rescale to unit diagonal).
nearest_correlation_matrix <- function(R) {
  eig <- eigen(R, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-8)
  M <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  d <- sqrt(diag(M))
  C <- M / outer(d, d)
  diag(C) <- 1
  C
}

is_valid_corr <- function(R) {
  ch <- tryCatch(chol(R), error = function(e) NULL)
  !is.null(ch)
}

# ---------- Build raw (possibly invalid) correlation matrix at time t ----------
build_raw_corr_matrix <- function(t, n_vars, pairs, pair_params) {
  R <- diag(n_vars)
  for (idx in seq_along(pairs)) {
    p <- pairs[[idx]]
    pp <- pair_params[[idx]]
    freq <- 2 * pi / pp$period
    r <- pp$base + pp$amp * sin(freq * t + pp$phase)
    r <- max(-0.95, min(0.95, r))
    R[p[1], p[2]] <- r
    R[p[2], p[1]] <- r
  }
  R
}

# ---------- Full dataset generator ----------
# Returns list(dataset = data.frame, corr_by_time = list of matrices,
#              n_corrections = integer)
generate_dataset <- function(n_subjects, n_time, n_vars, seed, noise_sd,
                             var_params, pair_params) {
  set.seed(seed)
  var_names <- var_names_for(n_vars)
  pairs <- pairs_for(n_vars)

  corr_by_time <- vector("list", n_time)
  chol_by_time <- vector("list", n_time)
  n_corrections <- 0

  for (t in 0:(n_time - 1)) {
    if (n_vars >= 2) {
      R <- build_raw_corr_matrix(t, n_vars, pairs, pair_params)
      if (!is_valid_corr(R)) {
        n_corrections <- n_corrections + 1
        R <- nearest_correlation_matrix(R)
      }
      L <- t(chol(R))  # lower-triangular Cholesky factor
    } else {
      R <- matrix(1, 1, 1)
      L <- matrix(1, 1, 1)
    }
    corr_by_time[[t + 1]] <- R
    chol_by_time[[t + 1]] <- L
  }

  rows <- vector("list", n_subjects * n_time)
  i <- 1
  for (subj in 1:n_subjects) {
    for (t in 0:(n_time - 1)) {
      z <- rnorm(n_vars)
      L <- chol_by_time[[t + 1]]
      correlated <- as.numeric(L %*% z)

      values <- numeric(n_vars)
      for (v in 1:n_vars) {
        vp <- var_params[[v]]
        trend  <- vp$intercept + vp$trendSlope * t
        season <- vp$seasonAmp * sin(2 * pi * t / vp$seasonPeriod)
        values[v] <- trend + season + correlated[v] * vp$varSd * noise_sd
      }

      row <- c(list(subject = subj, time = t), as.list(round(values, 4)))
      names(row)[3:(2 + n_vars)] <- var_names
      rows[[i]] <- row
      i <- i + 1
    }
  }

  dataset <- bind_rows(rows)
  list(dataset = dataset, corr_by_time = corr_by_time, n_corrections = n_corrections,
       n_vars = n_vars, var_names = var_names, pairs = pairs)
}

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background-color: #0f172a; color: #e2e8f0; }
    .well { background-color: #1e293b; border: 1px solid #334155; }
    h2, h3, h4 { color: #e2e8f0; }
    .form-control { background-color: #0f172a; color: #e2e8f0; border-color: #475569; }
    table.dataTable { color: #e2e8f0 !important; }
    .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter,
    .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_paginate {
      color: #e2e8f0 !important;
    }
  "))),

  titlePanel("Longitudinal Multivariate Dataset Generator"),
  p("Configurable number of correlated variables, time points, and subjects, with correlations drifting sinusoidally over time.",
    style = "color:#94a3b8;"),

  fluidRow(
    column(3,
           wellPanel(
             h4("Global settings"),
             numericInput("nVars", "# Variables", value = 5, min = 1, max = 15, step = 1),
             actionButton("applyDims", "Apply variable count", class = "btn-secondary", width = "100%"),
             p("Click this after changing # Variables to rebuild the parameter tables below.",
               style = "font-size:11px; color:#94a3b8; margin-top:6px;"),
             numericInput("nSubjects", "# Subjects", value = 20, min = 1, max = 5000),
             numericInput("nTime", "# Time points", value = 100, min = 2, max = 2000),
             numericInput("seed", "Random seed", value = 42),
             numericInput("noiseSd", "Noise scale", value = 1, min = 0.1, step = 0.1),
             hr(),
             actionButton("generate", "Generate dataset", class = "btn-primary", width = "100%"),
             br(), br(),
             downloadButton("downloadCsv", "Download CSV", width = "100%"),
             br(), br(),
             uiOutput("summaryNote"),
             uiOutput("correctionNote")
           )
    ),
    column(9,
           tabsetPanel(
             tabPanel("Variable parameters",
                      br(),
                      p("Per-variable trend, seasonality, and noise scale.", style="color:#94a3b8;"),
                      uiOutput("varParamInputs")
             ),
             tabPanel("Correlation parameters",
                      br(),
                      p("r(t) = base + amplitude x sin(2*pi*t/period + phase), clamped to [-0.95, 0.95].",
                        style="color:#94a3b8;"),
                      uiOutput("pairParamInputs")
             ),
             tabPanel("Correlation trajectories",
                      br(),
                      plotOutput("corrPlot", height = "450px")
             ),
             tabPanel("Sample subject trajectory",
                      br(),
                      numericInput("previewSubject", "Subject #", value = 1, min = 1, step = 1),
                      plotOutput("subjectPlot", height = "400px")
             ),
             tabPanel("Data preview",
                      br(),
                      DTOutput("dataTable")
             )
           )
    )
  )
)

# ============================================================
# Server
# ============================================================
server <- function(input, output, session) {

  # n_vars_active tracks the variable count actually used to build the
  # parameter UI. It only changes when "Apply variable count" is clicked,
  # so typing a new number doesn't destroy inputs mid-edit.
  n_vars_active <- reactiveVal(5)

  observeEvent(input$applyDims, {
    req(input$nVars)
    n_vars_active(max(1, min(15, round(input$nVars))))
  })

  # ---- Dynamic parameter inputs (rebuilt only when n_vars_active changes) ----
  output$varParamInputs <- renderUI({
    n_vars <- n_vars_active()
    var_names <- var_names_for(n_vars)
    tagList(lapply(1:n_vars, function(i) {
      default_slope <- 0.02 * i * ifelse(i %% 2 == 0, 1, -1)
      wellPanel(
        h4(var_names[i], style = paste0("color:", var_color_for(i, n_vars), ";")),
        fluidRow(
          column(2, numericInput(paste0("intercept_", i), "Intercept", value = 0, step = 0.1)),
          column(3, numericInput(paste0("trendSlope_", i), "Trend slope", value = default_slope, step = 0.01)),
          column(2, numericInput(paste0("seasonAmp_", i), "Season amplitude", value = 1 + (i - 1) * 0.3, step = 0.1)),
          column(3, numericInput(paste0("seasonPeriod_", i), "Season period", value = 20 + (i - 1) * 5, step = 1)),
          column(2, numericInput(paste0("varSd_", i), "SD scale", value = 1, min = 0.1, step = 0.1))
        )
      )
    }))
  })

  output$pairParamInputs <- renderUI({
    n_vars <- n_vars_active()
    var_names <- var_names_for(n_vars)
    pairs <- pairs_for(n_vars)
    if (length(pairs) == 0) {
      return(p("Only one variable selected - no pairwise correlations to configure.",
               style = "color:#94a3b8;"))
    }
    tagList(lapply(seq_along(pairs), function(idx) {
      p <- pairs[[idx]]
      default_base   <- 0.1 + 0.1 * ((idx - 1) %% 3)
      default_amp    <- 0.3 + 0.1 * ((idx - 1) %% 2)
      default_period <- 100 / (1 + ((idx - 1) %% 3))
      default_phase  <- (idx - 1) * 0.7
      wellPanel(
        h4(pair_label(p, var_names), style = "color:#cbd5e1;"),
        fluidRow(
          column(3, numericInput(paste0("base_", idx), "Base", value = default_base, min = -0.9, max = 0.9, step = 0.05)),
          column(3, numericInput(paste0("amp_", idx), "Amplitude", value = default_amp, min = 0, max = 0.95, step = 0.05)),
          column(3, numericInput(paste0("period_", idx), "Period", value = default_period, min = 2, step = 1)),
          column(3, numericInput(paste0("phase_", idx), "Phase", value = default_phase, step = 0.1))
        )
      )
    }))
  })

  # ---- Collect parameters reactively ----
  var_params <- reactive({
    n_vars <- n_vars_active()
    req(input[[paste0("intercept_", n_vars)]])  # wait until dynamic UI has rendered
    lapply(1:n_vars, function(i) {
      list(
        intercept    = input[[paste0("intercept_", i)]],
        trendSlope   = input[[paste0("trendSlope_", i)]],
        seasonAmp    = input[[paste0("seasonAmp_", i)]],
        seasonPeriod = input[[paste0("seasonPeriod_", i)]],
        varSd        = input[[paste0("varSd_", i)]]
      )
    })
  })

  pair_params <- reactive({
    n_vars <- n_vars_active()
    pairs <- pairs_for(n_vars)
    if (length(pairs) == 0) return(list())
    req(input[[paste0("base_", length(pairs))]])
    lapply(seq_along(pairs), function(idx) {
      list(
        base   = input[[paste0("base_", idx)]],
        amp    = input[[paste0("amp_", idx)]],
        period = input[[paste0("period_", idx)]],
        phase  = input[[paste0("phase_", idx)]]
      )
    })
  })

  # ---- Generate dataset (on button click) ----
  result <- eventReactive(input$generate, {
    n_vars <- n_vars_active()
    generate_dataset(
      n_subjects  = input$nSubjects,
      n_time      = input$nTime,
      n_vars      = n_vars,
      seed        = input$seed,
      noise_sd    = input$noiseSd,
      var_params  = var_params(),
      pair_params = pair_params()
    )
  }, ignoreNULL = FALSE)

  output$summaryNote <- renderUI({
    res <- tryCatch(result(), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    div(style = "color:#94a3b8; font-size:12px; margin-top:4px;",
        sprintf("Generated: %d variables x %d time points x %d subjects = %d rows.",
                res$n_vars, input$nTime, input$nSubjects, nrow(res$dataset)))
  })

  output$correctionNote <- renderUI({
    res <- tryCatch(result(), error = function(e) NULL)
    if (is.null(res) || res$n_corrections == 0) return(NULL)
    div(style = "color:#fbbf24; background:#78350f33; border:1px solid #92400e; border-radius:6px; padding:8px; font-size:12px; margin-top:6px;",
        sprintf("Note: %d of %d time points had inconsistent pairwise correlations and were automatically corrected to the nearest valid correlation matrix.",
                res$n_corrections, input$nTime))
  })

  output$corrPlot <- renderPlot({
    res <- result()
    req(res)
    if (length(res$pairs) == 0) {
      return(ggplot() + annotate("text", x = 0, y = 0, label = "Only one variable - no correlations to plot",
                                 color = "#94a3b8", size = 5) + theme_void())
    }
    var_names <- res$var_names
    corr_df <- bind_rows(lapply(seq_along(res$corr_by_time), function(t_idx) {
      R <- res$corr_by_time[[t_idx]]
      vals <- sapply(res$pairs, function(p) R[p[1], p[2]])
      setNames(as.data.frame(t(vals)), sapply(res$pairs, pair_label, var_names = var_names)) |>
        mutate(time = t_idx - 1)
    }))
    corr_long <- pivot_longer(corr_df, -time, names_to = "pair", values_to = "correlation")

    ggplot(corr_long, aes(x = time, y = correlation, color = pair)) +
      geom_line(linewidth = 0.8) +
      ylim(-1, 1) +
      labs(x = "Time", y = "Correlation", color = "Pair") +
      theme_minimal(base_size = 13) +
      theme(
        plot.background = element_rect(fill = "#1e293b", color = NA),
        panel.background = element_rect(fill = "#1e293b", color = NA),
        panel.grid = element_line(color = "#334155"),
        text = element_text(color = "#e2e8f0"),
        axis.text = element_text(color = "#94a3b8"),
        legend.background = element_rect(fill = "#1e293b"),
        legend.text = element_text(color = "#e2e8f0")
      )
  })

  output$subjectPlot <- renderPlot({
    res <- result()
    req(res)
    subj_id <- max(1, min(input$previewSubject, input$nSubjects))
    df <- res$dataset |> filter(subject == subj_id)
    df_long <- pivot_longer(df, all_of(res$var_names), names_to = "variable", values_to = "value")
    color_map <- setNames(sapply(1:res$n_vars, var_color_for, n_vars = res$n_vars), res$var_names)

    ggplot(df_long, aes(x = time, y = value, color = variable)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = color_map) +
      labs(x = "Time", y = "Value", color = "Variable",
           title = paste("Subject", subj_id)) +
      theme_minimal(base_size = 13) +
      theme(
        plot.background = element_rect(fill = "#1e293b", color = NA),
        panel.background = element_rect(fill = "#1e293b", color = NA),
        panel.grid = element_line(color = "#334155"),
        text = element_text(color = "#e2e8f0"),
        axis.text = element_text(color = "#94a3b8"),
        legend.background = element_rect(fill = "#1e293b"),
        legend.text = element_text(color = "#e2e8f0")
      )
  })

  output$dataTable <- renderDT({
    res <- result()
    req(res)
    datatable(res$dataset, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$downloadCsv <- downloadHandler(
    filename = function() {
      res <- tryCatch(result(), error = function(e) NULL)
      if (is.null(res)) return("longitudinal_dataset.csv")
      sprintf("longitudinal_dataset_%dvars_%dtime_%dsubj.csv",
              res$n_vars, input$nTime, input$nSubjects)
    },
    content = function(file) {
      res <- result()
      req(res)
      write.csv(res$dataset, file, row.names = FALSE)
    }
  )
}

# ============================================================
# Run app
# ============================================================
shinyApp(ui = ui, server = server)
