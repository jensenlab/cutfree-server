#!/usr/local/bin/Rscript

library(magrittr)

library(loggr)
library(cutfree)


# ========= helper functions =============
`%<in>%` <- function(x,y) all(x %in% y) & all(y %in% x)

# converts oligos of the form
#     "AGT{G,C}G{T,A}"
# to the corresponding IUB-coded string
#     "AGTSGW"
sets_to_codes <- function(string, codes=cutfree::IUB_CODES) {
  blocks <- (string %>%
    stringr::str_split("[\\{\\}]"))[[1]] %>% 
    Filter(function(x) nchar(x) > 0, .)
  for (i in seq_along(blocks)) {
    if (stringr::str_detect(blocks[i], ",")) {
      bases <- stringr::str_split(blocks[i], "[^ACGT]")[[1]] %>% 
        Filter(function(x) nchar(x) > 0, .)
      blocks[i] <- cutfree::which_name(
        sapply(codes, `%<in>%`, bases)
      )
    }
  }
  paste0(blocks, collapse="")
}

codes_to_sets <- function(string, codes=cutfree::IUB_CODES) {
  to_set <- function(x) {
    if (length(x) == 1) x else paste0("{", paste(x, collapse=","), "}")
  }
  cutfree::char2str(lapply(codes, to_set)[cutfree::str_to_vector(string)])
}

# ========= main functions =============

say_error <- function(title) {
  sayf("*** ERROR: %s ***\n", title)
  indent_log()
}

create_safe_randomize_oligo <- function(codes=IUB_CODES,
                                        max_oligo_length=30,
                                        max_blocks=2,
                                        max_time=300,
                                        max_sites=3) {

  safe_randomize_oligo <- function(unsafe_starting_oligo,
                                   unsafe_sites,
                                   unsafe_min_blocks=1,
                                   unsafe_re_randomize=FALSE,
                                   unsafe_obj_frac=1.0,
                                   raw=FALSE) {
    in_codes <- function(strs, cds=codes) {
      aux <- function(s) all(stringr::str_split(s, "")[[1]] %in% names(cds))
      all(sapply(strs, aux))
    }

    # convert CGI parameters
    unsafe_starting_oligo %<>% toupper()
    use_set_format <-  nchar(sets_to_codes(unsafe_starting_oligo)) != nchar(unsafe_starting_oligo)
    unsafe_starting_oligo %<>% sets_to_codes()
    unsafe_obj_frac %<>% as.numeric()
    unsafe_sites <- stringr::str_split(httpuv::decodeURI(unsafe_sites), ",\\s*")[[1]] %>%
      toupper()
    unsafe_re_randomize %<>% as.logical()
    unsafe_min_blocks %<>% as.integer()
    raw %<>% as.logical()

    safe_randomize_oligo_aux <- function() {
      if (!raw)
        sayf("<pre>\n")

      if (nchar(unsafe_starting_oligo) > max_oligo_length) {
        say_error("Oligo too long")
        sayf("The CutFree server allows a maximum oligo length of %i bases.", max_oligo_length)
        sayf("The oligo you provided is %i bases.\n", nchar(unsafe_starting_oligo))
        say("To use CutFree on longer oligos, please download and install the software locally.")
        unindent_log()
        return()
      }

      if (!in_codes(unsafe_starting_oligo)) {
        say_error("Unrecognized Code")
        say("The starting oligo contains a letter that is not a valid IUB code.")
        unindent_log()
        return()
      }

      if (!in_codes(unsafe_sites)) {
        say_error("Unrecognized Code")
        sayf("The restriction sites contain a letter that is not a valid IUB code.")
        unindent_log()
        return()
      }

      if (unsafe_min_blocks > max_blocks) {
        say_error("Too many blocks requested")
        sayf("The CutFree server allows a maximum of %i blocks.", max_blocks)
        sayf("You requested %i blocks.\n", unsafe_min_blocks)
        say("To use CutFree with more blocks, please download and install the software locally.")
        unindent_log()
        return()
      }

      if (length(unsafe_sites) > max_sites) {
        say_error("Too many sites requested")
        sayf("The CutFree server allows a maximum of %i sites", max_sites)
        sayf("You requested %i sites.\n", length(unsafe_sites))
        say("To use CutFree with more sites, please download and install the software locally.")
        unindent_log()
        return()
      }

      result <- cutfree(starting_oligo=unsafe_starting_oligo,
                        min_blocks=unsafe_min_blocks,
                        sites=unsafe_sites,
                        re_randomize=unsafe_re_randomize,
                        obj_frac=unsafe_obj_frac,
                        quiet=TRUE,
                        maxtime=max_time)

      if (use_set_format)
        sayf("Input Oligo: \n\n   %s\n", unsafe_starting_oligo %>% codes_to_sets())
      else
        sayf("Input Oligo: \n\n   %s\n", unsafe_starting_oligo)
      
      if (use_set_format)
        sayf("Blocked Oligo:\n\n   %s\n", result$code %>% codes_to_sets())
      else
        sayf("Blocked Oligo:\n\n   %s\n", result$code)
      
      say("Base Distribution:\n")
      print_oligo_block(result$code, indent="   ")

      sayf("\n\nSolution status: %s.", result$sol_orig$status)
      sayf("Solver used %i iterations in %g seconds.", result$sol_orig$itercount, result$sol_orig$runtime)
      sayf("Final objective value was %g.", result$sol_orig$objval)
    }

    return(capture.output(safe_randomize_oligo_aux()) %>% paste(collapse="\n"))
  }

  return(safe_randomize_oligo)
}


start_oligo_server <- function(randomizer=create_safe_randomize_oligo(),
                               ip="127.0.0.1") {
  jug::jug() %>%
    jug::get(path="/cutfree", jug::decorate(randomizer)) %>%
    jug::simple_error_handler() %>%
    jug::serve_it(host=ip)
}

start_oligo_server(ip="128.174.125.122")




