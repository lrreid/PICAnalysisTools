''''
Functions for converting units between any two orders of magnitude

'''

def order_of_magnitude(unit):
     match unit:
          case "quecto" | "q":
               return 1e-30
          case "ronto" | "r":
               return 1e-27
          case "yocto" | "y":
               return 1e-24
          case "zepto" | "z":
               return 1e-21
          case "atto" | "a":
               return 1e-18
          case "femto" | "f":
               return 1e-15
          case "pico" | "p":
               return 1e-12
          case "nano" | "n":
               return 1e-9
          case "micro" | "u":
               return 1e-6
          case "milli" | "m":
               return 1e-3
          case "centi" | "c":
               return 1e-2
          case "one" | "":
               return 1
          case "hecto" | "h":
               return 1e2
          case "kilo" | "k":
               return 1e3
          case "mega" | "Mega" | "M":
               return 1e6
          case "giga" | "Giga" | "G":
               return 1e9
          case "tera" | "Tera" | "T":
               return 1e12
          case "peta" | "Peta" | "P":
               return 1e15
          case "exa" | "Exa" | "E":
               return 1e18
          case "zetta" | "Zetta" | "Z":
               return 1e21
          case "yotta" | "Yotta" | "Y":
               return 1e24
          case "ronna" | "Ronna" | "R":
               return 1e27
          case "quetta" | "Quetta" | "Q":
               return 1e30


def magnitude_conversion(values, unit_1, unit_2, reciprocal_units = False):

     if unit_1 == unit_2:
          return values
     else:
          order_1 = order_of_magnitude(unit_1)
          order_2 = order_of_magnitude(unit_2)

          if reciprocal_units == True:
               values_new = (order_2/order_1) * values
          elif reciprocal_units == False:
               values_new = (order_1/order_2) * values

          return values_new

def magnitude_conversion_area(values, unit_1, unit_2, reciprocal_units = False):

     if unit_1 == unit_2:
          return values
     else:
          order_1 = order_of_magnitude(unit_1)
          order_2 = order_of_magnitude(unit_2)

          if reciprocal_units == True:
               values_new = (order_2/order_1)**2 * values
          elif reciprocal_units == False:
               values_new = (order_1/order_2)**2 * values

          return values_new     

def magnitude_conversion_vol(values, unit_1, unit_2, reciprocal_units = False):

     if unit_1 == unit_2:
          return values
     else:
          order_1 = order_of_magnitude(unit_1)
          order_2 = order_of_magnitude(unit_2)

          if reciprocal_units == True:
               values_new = (order_2/order_1)**3 * values
          elif reciprocal_units == False:
               values_new = (order_1/order_2)**3 * values

          return values_new
     

def get_order_letter(unit):
     match unit:
          case "quecto" | "q":
               return "q"
          case "ronto" | "r":
               return "r"
          case "yocto" | "y":
               return "y"
          case "zepto" | "z":
               return "z"
          case "atto" | "a":
               return "a"
          case "femto" | "f":
               return "f"
          case "pico" | "p":
               return "p"
          case "nano" | "n":
               return "n"
          case "micro" | "u":
               return "u"
          case "milli" | "m":
               return "m"
          case "centi" | "c":
               return "c"
          case "one" | "":
               return ""
          case "hecto" | "h":
               return "h"
          case "kilo" | "k":
               return "k"
          case "mega" | "Mega" | "M":
               return "M"
          case "giga" | "Giga" | "G":
               return "G"
          case "tera" | "Tera" | "T":
               return "T"
          case "peta" | "Peta" | "P":
               return "P"
          case "exa" | "Exa" | "E":
               return "E"
          case "zetta" | "Zetta" | "Z":
               return "Z"
          case "yotta" | "Yotta" | "Y":
               return "Y"
          case "ronna" | "Ronna" | "R":
               return "R"
          case "quetta" | "Quetta" | "Q":
               return "Q"