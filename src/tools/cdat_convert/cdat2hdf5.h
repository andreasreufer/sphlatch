#include <string>

std::pair<std::string, std::string> FloatOptParser(const std::string &str)
{
  if (str == "--float")
    {
      return std::make_pair(std::string("float"), std::string("true"));
    }
  else
    {
      return std::make_pair(std::string(), std::string());
    }
}

std::pair<std::string, std::string> DoubleOptParser(const std::string &str)
{
  if (str == "--double")
    {
      return std::make_pair(std::string("double"), std::string("true"));
    }
  else
    {
      return std::make_pair(std::string(), std::string());
    }
}

// The last extra parser never seems to work, so insert a dummy:

std::pair<std::string, std::string> DummyOptParser(const std::string &str)
{
  return std::make_pair(std::string(), std::string());
}

