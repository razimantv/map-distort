#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

const double distortion_threshold = 0.01, outer_threshold = 2,
             accuracy_threshold = 1.01;
const double pi = acos(-1);

// A pair of (x,y) coordinates is a point
typedef std::pair<double, double> point;
// A simply connected non-intersecting polygon is a piece
typedef std::vector<point> piece;   
// A group of non-intersecting pieces is a region (state/district...)
typedef std::vector<piece> region;
// And a group of regions forms a map
typedef std::vector<region> map;

// Find area and centroid of piece from line integral using Green's theorem
std::pair<double, point> area_centroid(const piece& p) {
  int N = p.size();
  double intdA = 0, intxdA = 0, intydA = 0;
  for (int i = 1; i < N; ++i) {
    auto [xp, yp] = p[i - 1];
    auto [x, y] = p[i];
    intdA += y * xp - x * yp;
    intxdA += (y - yp) * (x * x + x * xp + xp * xp) / 3;
    intydA += (xp - x) * (y * y + y * yp + yp * yp) / 3;
  }
  return {std::abs(intdA) / 2, {intxdA / intdA, intydA / intdA}};
}

// Now propagate it to regions and maps
template <class T>
std::pair<double, point> area_centroid(const T& r) {
  double intdA = 0, intxdA = 0, intydA = 0;
  for (auto& p : r) {
    if(p.empty()) continue;
    auto [dA, dxy] = area_centroid(p);
    intdA += dA;
    intxdA += dA * dxy.first;
    intydA += dA * dxy.second;
  }
  return {intdA, {intxdA / intdA, intydA / intdA}};
}

// Find the convex hull of a region using Graham's scan
piece convexhull(const region& r) {
  const double eps = 1e-9;
  point P0 = r[0][0];
  auto& [x0, y0] = P0;
  for (const piece& p : r)
    for (auto [x, y] : p) {
      if (y < y0 or (y == y0 and x < x0)) x0 = x, y0 = y;
    }

  piece pvec, hull;
  for (const piece& p : r)
    for (auto [x, y] : p) {
      if (std::abs(x - x0) > eps or std::abs(y - y0) > eps)
        pvec.push_back({x, y});
    }

  auto comp = [eps](const point& p0, const point& p1, const point& p2) {
    return (p1.first - p0.first) * (p2.second - p0.second) -
               (p2.first - p0.first) * (p1.second - p0.second) <
           -eps;
  };
  sort(pvec.begin(), pvec.end(),
       [&P0, &comp](const point& p1, const point& p2) {
         return comp(P0, p1, p2);
       });

  hull.push_back(P0);
  pvec.push_back(P0);
  int H = 1;
  for (auto& p : pvec) {
    while (H > 1 and !comp(hull[H - 2], hull[H - 1], p)) hull.pop_back(), --H;
    hull.push_back(p);
    ++H;
  }
  return hull;
}

// Poor implementation of a subset of svg path reading
// Only considers piecewise linear paths
map readsvgpaths(const std::string& filename) {
  std::ifstream in(filename);
  map ret;
  for (std::string s; std::getline(in, s);) {
    ret.push_back({});
    auto& subreg = ret.back();
    std::istringstream iss(s);

    double x = 0, y = 0, x0, y0;
    for (char motion; iss;) {
      if (!(iss >> std::ws)) break;
      int next = iss.peek();
      if (next == EOF)
        break;
      else if (isalpha(next))
        iss >> motion;
      std::cerr << motion;

      if (motion == 'm' or motion == 'l') {
        double dx, dy;
        iss >> dx >> dy;
        x += dx;
        y += dy;
        if (motion == 'm') {
          x0 = x;
          y0 = y;
          motion = 'l';
          subreg.push_back({});
        }
      } else if (motion == 'M' or motion == 'L') {
        iss >> x >> y;
        if (motion == 'M') {
          x0 = x;
          y0 = y;
          motion = 'L';
          subreg.push_back({});
        }
      } else if (motion == 'h') {
        double dx;
        iss >> dx;
        x += dx;
      } else if (motion == 'H') {
        iss >> x;
      } else if (motion == 'v') {
        double dy;
        iss >> dy;
        y += dy;
      } else if (motion == 'V') {
        iss >> y;
      } else if (motion == 'z') {
        x = x0;
        y = y0;
      } else {
        std::cerr << "Error " << motion << "\n";
        exit(1);
      }
      subreg.back().push_back({x, y});
    }
    std::cerr << "\n";
  }
  return ret;
}

// Read and normalise the weights for each part of the map
std::vector<double> readweights(const std::string& filename) {
  std::ifstream in(filename);
  std::vector<double> ret;
  double p, tot = 0;
  while (in >> p) {
    ret.push_back(p);
    tot += p;
  }
  for (double& p : ret) p /= tot;
  return ret;
}

// Aspect ratio of the map
double aspect(const map& M) {
  double xmin = std::numeric_limits<double>::max(),
         xmax = std::numeric_limits<double>::min(), ymin = xmin, ymax = xmax;
  for (auto& reg : M)
    for (auto& p : reg)
      for (auto& [x, y] : p) {
        xmax = std::max(x, xmax);
        xmin = std::min(x, xmin);
        ymax = std::max(y, ymax);
        ymin = std::min(y, ymin);
      }

  return (ymax - ymin) / (xmax - xmin);
}

// Move the origin to the centroid of the map
void recentre(map& M) {
  auto [A, centre] = area_centroid(M);
  auto [cx, cy] = centre;
  for (auto& reg : M)
    for (auto& p : reg)
      for (auto& [x, y] : p) x -= cx, y -= cy;
}

// Major work: Stretch or compress a region by a factor
// without messing up rest of the map
void distort(map& M, int id, double distortion_factor, double aspect_target) {
  // Start with the centroid of the convex hull
  auto H = convexhull(M[id]);
  auto [A, centre] = area_centroid(H);
  auto [cx, cy] = centre;

  // Sort the points on the convex hull using angle around centroid
  // Unnecessary log N factor here due to map, TODO: Optimise if too slow
  std::map<double, point> anglehull;
  for (auto [x, y] : H) {
    double theta = atan2(y - cy, x - cx);
    anglehull[theta] = {x - cx, y - cy};
  }

  // Close the "gap" in the hull near -pi/pi
  {
    auto [theta1, c1] = *anglehull.begin();
    auto [theta2, c2] = *anglehull.rbegin();
    anglehull[theta1 + 2 * pi] = c1;
    anglehull[theta2 - 2 * pi] = c2;
  }

  // Express every line as a*x+b=c, such that c=1 for the convex hull
  std::vector<std::tuple<double, double, double>> theta_ab;
  for (auto mit = anglehull.begin();;) {
    auto [x1, y1] = mit->second;
    double theta1 = mit->first;
    if(++mit == anglehull.end()) break;
    auto [x2, y2] = mit->second;
    double theta2 = mit->first;
    if(theta2 - theta1 < 1e-6) continue;
    double det = x1 * y2 - x2 * y1, a = (y2-y1) / det, b = (x1-x2) / det;
    theta_ab.push_back({theta2, a, b});
    //std::cerr << theta2 << " " << a << " " << b << "\n";
  }

  // Stretch each point in the map depending on the angular region
  // Mapping c in [0,1] to [0, distortion_factor]
  //          and [1, outer_threshold] to [distortion_factor, outer_threshold]
  for (auto& reg : M)
    for (auto& p : reg)
      for (auto& [x, y] : p) {
        double dx = x - cx, dy = y - cy;
        double theta = atan2(dy, dx);

        // Find angular sector the point falls in
        auto [thetahull, a, b] = *upper_bound(theta_ab.begin(), theta_ab.end(),
                                              std::make_tuple(theta, 0., 0.));

        // c value range tells us how to distort
        double c = a * dx + b * dy;
        if(c > outer_threshold) continue;

        double cnew = 1;
        if(c < 1) {
          cnew = c * distortion_factor;
        } else {
          double m1 = (outer_threshold - distortion_factor) /
                      (outer_threshold - 1),
                 m2 = outer_threshold * (distortion_factor - 1) /
                      (outer_threshold - 1);
          cnew = m1 * c + m2;
        }
        x = cx + dx * cnew / c;
        y = cy + dy * cnew / c;
        //std::cerr << dx << " "<< dy << " " << c << " " << cnew << "\n";
      }

  // Points have been stretched, now recentre the map
  recentre(M);

  // And fix its aspect ratio without affecting the area
  double aspect_modifier = sqrt(aspect(M) / aspect_target);
  for (auto& reg : M)
    for (auto& p : reg)
      for (auto& [x, y] : p) x *= aspect_modifier, y /= aspect_modifier;
}

int main() {
  std::string mapfilename, popfilename;
  int rounds;
  std::cin >> mapfilename >> popfilename >> rounds;
  
  map M = readsvgpaths(mapfilename);
  auto weights = readweights(popfilename);
  int N = M.size();

  for (int i = 0; i < N; ++i)
    if (weights[i] == 0) M[i].clear();

  recentre(M);
  double aspect_target = aspect(M);

  std::vector<int> cnt(N);
  for (int i = 0; ; ++i) {
    // Normalise total area of the map to 1
    std::vector<double> areas;
    std::transform(M.begin(), M.end(), std::back_inserter(areas),
                   [](const region& d) { return area_centroid(d).first; });

    double totarea = std::accumulate(areas.begin(), areas.end(), 0.);
    double rescale_factor = sqrt(totarea);
    for (auto& reg : M)
      for (auto& subreg : reg)
        for (auto& xy : subreg) xy = {xy.first / rescale_factor, xy.second / rescale_factor};
    for (double& a : areas) a /= totarea;

    // Find the most deviating region, but reduce its weight
    // if we have processed it too many times before
    double worst = 0, absworst = 0;
    int worst_id = -1;
    for (int i = 0; i < N; ++i) {
      if(weights[i] == 0) continue;
      double bad = areas[i] / weights[i];
      bad = std::max(bad, 1 / bad);
      absworst = std::max(absworst, bad);
      bad = (bad - 1) / (cnt[i] + 1);
      if (bad > worst) {
        worst = bad;
        worst_id = i;
      }
    }
    // Find the amount of distortion needed to "correct" the area
    double distortion_factor = sqrt(weights[worst_id] / areas[worst_id]);

    // But pin it to threshold level
    if (distortion_factor > 1 + distortion_threshold)
      distortion_factor = 1 + distortion_threshold;
    else if (distortion_factor < 1 - distortion_threshold)
      distortion_factor = 1 - distortion_threshold;

    std::cerr << worst_id << " " << cnt[worst_id] << " " << areas[worst_id]
              << " " << weights[worst_id] << " " << absworst << " "
              << distortion_factor << "\n";

    // Stop if all regions are fixed to desired accuracy
    if (absworst < accuracy_threshold) break;

    if (i == rounds) break;
    // Distort the map
    distort(M, worst_id, distortion_factor, aspect_target);
    ++cnt[worst_id];
  }

  // Output the map with region ID in gnuplot-readable format
  for (int i = 0; i < N; ++i) {
    for (auto& r : M[i]) {
      for (auto [x, y] : r) std::cout << x << " " << y << " " << i + 1 << "\n";
      std::cout << "\n";
    }
  }

  return 0;
}
