library(shiny)
library(ssadequacy)

ui <- fluidPage(
  titlePanel("Sample Size Adequacy (SSA) planner"),
  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("type", "Design", c("Confirmatory factor analysis" = "cfa",
                                      "Multiple regression" = "regression", "One-way ANOVA" = "anova")),
      conditionalPanel("input.type == 'cfa'",
        numericInput("nf", "Factors", 3, 1, 10), numericInput("nind", "Indicators per factor", 5, 3, 12),
        sliderInput("lam", "Standardized loadings (range)", .2, .9, c(.55, .60), .01),
        numericInput("phi", "Factor correlation", .30, -.9, .9, .05)),
      conditionalPanel("input.type == 'regression'",
        numericInput("k", "Predictors", 3, 1, 30), numericInput("R2", "Population R-squared", .10, .01, .9, .01),
        numericInput("rho", "Predictor intercorrelation", .30, 0, .9, .05)),
      conditionalPanel("input.type == 'anova'",
        numericInput("g", "Groups", 2, 2, 10), numericInput("f", "Cohen's f", .25, .02, 1, .01)),
      checkboxInput("nn", "Non-normal data (skewness 2, excess kurtosis 7; MLR for CFA)", FALSE),
      hr(),
      numericInput("delta", "Target 95% CI width (delta)", .20, .01, 2, .01),
      sliderInput("power", "Target power", .5, .99, .80, .01),
      numericInput("alpha", "Alpha (two-sided)", .05, .001, .2, .005),
      checkboxInput("use_ass", "Require assurance for the width", FALSE),
      conditionalPanel("input.use_ass", sliderInput("ass", "Assurance", .5, .99, .80, .01)),
      actionButton("go", "Analytic bottleneck analysis", class = "btn-primary"),
      hr(),
      numericInput("reps", "Monte Carlo replications per N", 300, 50, 5000, 50),
      actionButton("sim", "Verify by simulation (slower)")
    ),
    mainPanel(width = 8,
      h4("Analytic result"), verbatimTextOutput("analytic"),
      plotOutput("crossover", height = 330),
      h4("Justification text"), verbatimTextOutput("report"),
      h4("Simulation"), verbatimTextOutput("simout")
    )
  )
)

server <- function(input, output, session) {
  design <- reactive({
    sk <- if (input$nn) 2 else 0; ku <- if (input$nn) 7 else 0
    switch(input$type,
      cfa = ssa_design_cfa(input$nf, input$nind, loadings = input$lam, phi = input$phi, skew = sk, exkurt = ku),
      regression = ssa_design_regression(input$k, input$R2, input$rho, skew = sk, exkurt = ku),
      anova = ssa_design_anova(input$g, input$f, skew = sk, exkurt = ku))
  })
  res <- eventReactive(input$go, {
    d <- design(); v <- ssa_unit_variance(d, M = 5e4)
    list(d = d, v = v, a = ssa_analytic(d, input$delta, input$power, input$alpha,
                                        if (input$use_ass) input$ass else NULL, v = v))
  })
  output$analytic <- renderPrint(print(res()$a))
  output$report <- renderPrint(invisible(ssa_report(res()$a)))
  output$crossover <- renderPlot({
    r <- res(); th <- mean(abs(r$d$theta)); ds <- seq(.25, 3, length.out = 120) * th
    Np <- r$a$N["Detectability"]
    Nd <- vapply(ds, function(dl) ssa_analytic(r$d, dl, input$power, input$alpha, v = r$v)$N["Precision"], 0)
    plot(ds, Nd, log = "y", type = "l", lwd = 2, lty = 2, xlab = expression("Target CI width " * delta),
         ylab = "Required N (log scale)", ylim = range(c(Nd, Np)))
    abline(h = Np, lwd = 2); abline(v = r$a$r_star * th, lty = 3)
    lines(ds, pmax(Nd, Np), col = "#2d4f7c", lwd = 4)
    points(input$delta, r$a$N_star, pch = 19, cex = 1.4, col = "#c8742c")
    legend("topright", c("Detectability", "Precision", "N* (max)", "crossover r* x |theta|"),
           lty = c(1, 2, 1, 3), lwd = c(2, 2, 4, 1), col = c(1, 1, "#2d4f7c", 1), bty = "n")
  })
  output$simout <- renderPrint({
    input$sim
    isolate({
      if (input$sim == 0) return(cat("Press 'Verify by simulation' to run a Monte Carlo check."))
      s <- ssa_simulate(design(), input$delta, input$power, input$alpha,
                        if (input$use_ass) input$ass else NULL, reps = input$reps)
      print(s)
    })
  })
}

shinyApp(ui, server)
