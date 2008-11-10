#include <iostream>
#include <iomanip>
#include <sstream>

int main()
{
  std::cout << "#!/bin/bash\n"; 

  for (size_t i = 0; i < 897; i++)
  {
    std::string outStr;

    std::ostringstream idxSS;
    idxSS << i;

    for (size_t j = idxSS.str().size(); j < 3; j++)
      outStr.append("0");
    outStr.append( idxSS.str() );

    outStr.append(" ");

    const double time = (static_cast<int>(i) - 27)/36.;
    std::ostringstream timeSS;
    timeSS << std::fixed << std::setprecision(2) << time;
    
    for (size_t j = timeSS.str().size(); j < 5; j++)
      outStr.append(" ");
    outStr.append(timeSS.str());

    std::cout << "./annotate_frame.sh " << outStr << "\n";
  }
  std::cout << "gifsicle --delay 5 --colors 256 out*.gif   > anim.gif\n";
}
