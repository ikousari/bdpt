#ifndef BDPTCAMERA_H
#define BDPTCAMERA_H

//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
// Extended by Isaac Kousari in 2025 for Bidirectional Path Tracing

#include "hittable.h"
#include "pdf.h"
#include "material.h"

struct Vertex {
    Vertex(){}
    hit_record rec;

    vec3 wi;
    vec3 wo;
    bool isDelta = false;
    // pdf of sampling the *next* direction from this vertex, measured per solid angle
    double pdfFwd = 0.0;
    // pdf of sampling the *previous* direction (reverse) per solid angle
    double pdfRev = 0.0;
    // Cumulative throughput up to (and excluding) this vertex. Follows typical path tracer convention.
    color throughput = color(1,1,1);
};

struct Path {
    std::vector<Vertex> vertices;
    color throughput = color(1,1,1);
};

class bdptcamera {
  public:
    double aspect_ratio      = 1.0;  // Ratio of image width over height
    int    image_width       = 100;  // Rendered image width in pixel count
    int    samples_per_pixel = 10;   // Count of random samples for each pixel
    int    max_depth         = 10;   // Maximum number of ray bounces into scene
    color  background;               // Scene background color

    double vfov     = 90;              // Vertical view angle (field of view)
    point3 lookfrom = point3(0,0,0);   // Point camera is looking from
    point3 lookat   = point3(0,0,-1);  // Point camera is looking at
    vec3   vup      = vec3(0,1,0);     // Camera-relative "up" direction

    double defocus_angle = 0;  // Variation angle of rays through each pixel
    double focus_dist = 10;    // Distance from camera lookfrom point to plane of perfect focus

    void render(const hittable& world, const hittable& lights) {
        initialize();

        std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        for (int j = 0; j < image_height; j++) {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; i++) {
                color pixel_color(0,0,0);
                for (int s_j = 0; s_j < sqrt_spp; s_j++) {
                    for (int s_i = 0; s_i < sqrt_spp; s_i++) {
                        Path cameraPath = generate_camera_subpath(world,max_depth,i,j,s_i,s_j);
                        Path lightPath = generate_light_subpath(world,lights,max_depth);
                        pixel_color += ConnectPaths(world,lights,cameraPath,lightPath);
                    }
                }
                write_color(std::cout, pixel_samples_scale * pixel_color);
            }
        }

        std::clog << "\rDone.                 \n";
    }

  private:
    int    image_height;         // Rendered image height
    double pixel_samples_scale;  // Color scale factor for a sum of pixel samples
    int    sqrt_spp;             // Square root of number of samples per pixel
    double recip_sqrt_spp;       // 1 / sqrt_spp
    double viewport_height;
    double viewport_width;
    point3 center;               // Camera center
    point3 pixel00_loc;          // Location of pixel 0, 0
    vec3   pixel_delta_u;        // Offset to pixel to the right
    vec3   pixel_delta_v;        // Offset to pixel below
    vec3   u, v, w;              // Camera frame basis vectors
    vec3   defocus_disk_u;       // Defocus disk horizontal radius
    vec3   defocus_disk_v;       // Defocus disk vertical radius

    void initialize() {
        image_height = int(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        sqrt_spp = int(std::sqrt(samples_per_pixel));
        pixel_samples_scale = 1.0 / (sqrt_spp * sqrt_spp);
        recip_sqrt_spp = 1.0 / sqrt_spp;

        center = lookfrom;

        // Determine viewport dimensions.
        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta/2);
        viewport_height = 2 * h * focus_dist;
        viewport_width = viewport_height * (double(image_width)/image_height);

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        vec3 viewport_u = viewport_width * u;    // Vector across viewport horizontal edge
        vec3 viewport_v = viewport_height * -v;  // Vector down viewport vertical edge

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u/2 - viewport_v/2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        // Calculate the camera defocus disk basis vectors.
        auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    ray get_ray(int i, int j, int s_i, int s_j) const {
        // Construct a camera ray originating from the defocus disk and directed at a randomly
        // sampled point around the pixel location i, j for stratified sample square s_i, s_j.

        auto offset = sample_square_stratified(s_i, s_j);
        auto pixel_sample = pixel00_loc
                          + ((i + offset.x()) * pixel_delta_u)
                          + ((j + offset.y()) * pixel_delta_v);

        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction = pixel_sample - ray_origin;
        auto ray_time = random_double();

        return ray(ray_origin, ray_direction, ray_time);
    }

    vec3 sample_square_stratified(int s_i, int s_j) const {
        // Returns the vector to a random point in the square sub-pixel specified by grid
        // indices s_i and s_j, for an idealized unit square pixel [-.5,-.5] to [+.5,+.5].

        auto px = ((s_i + random_double()) * recip_sqrt_spp) - 0.5;
        auto py = ((s_j + random_double()) * recip_sqrt_spp) - 0.5;

        return vec3(px, py, 0);
    }

    vec3 sample_square() const {
        // Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
        return vec3(random_double() - 0.5, random_double() - 0.5, 0);
    }

    vec3 sample_disk(double radius) const {
        // Returns a random point in the unit (radius 0.5) disk centered at the origin.
        return radius * random_in_unit_disk();
    }

    point3 defocus_disk_sample() const {
        // Returns a random point in the camera defocus disk.
        auto p = random_in_unit_disk();
        return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
    }

    double pdf() const 
    {
        // Pixel area
        double pixelArea = (double)(image_width*image_height);
        return 1.0 / pixelArea;
    }

    double pdf(vec3 rayDir) const
    {
        // Camera space z-component (cosine of angle from view direction)
        double cosTheta = dot(unit_vector(rayDir), -w);
        
        if (cosTheta < 0.01)
            return 0.0;
        
        double num_pixels = (double)(image_width * image_height);
        double viewport_area = viewport_width * viewport_height;
        double cos_cubed = cosTheta * cosTheta * cosTheta;
        
        return num_pixels / (viewport_area * cos_cubed);
    }

    color ray_color(const ray& r, int depth, const hittable& world, const hittable& lights)
    const {
        return color(0,0,0);
    }

    Path generate_camera_subpath(const hittable &world, int maxDepth, int px, int py, int sx, int sy)
    {
        Path path;
        ray ray_i = get_ray(px, py, sx, sy);
        color throughput = color(1, 1, 1);
        double cameraPdf = pdf(unit_vector(ray_i.direction()));
        
        for (int depth = 0; depth < maxDepth; ++depth)
        {
            Vertex vertex;
            
            // Intersect ray with scene
            if (!world.hit(ray_i, interval(0.001, std::numeric_limits<double>::infinity()), vertex.rec))
                break;
            
            vertex.wi = unit_vector(-ray_i.direction());
            
            // Russian roulette
            if (depth > 3)
            {
                double continueProbability = fmin(0.95, color_to_luminance(throughput));
                if (random_double() > continueProbability)
                    break;
                throughput /= continueProbability;
            }

            vertex.throughput = throughput;
            BSDFSample bsdfSample = vertex.rec.mat->Sample(
                vertex.wi,
                vertex.rec,
                random_double(),
                random_double(),
                random_double()
            );
            
            if (bsdfSample.pdf == 0)
                break;
            
            // Store sampled direction and PDF
            vertex.wo = bsdfSample.dir;
            vertex.isDelta = bsdfSample.isDelta;
            vertex.pdfFwd = bsdfSample.pdf; 

            if (depth == 0)
            {
                vertex.pdfRev = cameraPdf;  // Camera sampling PDF
            }
            else
            {
                vertex.pdfRev = vertex.rec.mat->Pdf(vertex.wo, vertex.wi, vertex.rec);
            }
            
            path.vertices.push_back(vertex);
            
            if (false)
            {
                throughput = throughput * bsdfSample.val / bsdfSample.pdf;
            }
            else
            {
                double cosTheta = abs(dot(vertex.wo, vertex.rec.normal));
                throughput = throughput * bsdfSample.val * cosTheta / bsdfSample.pdf;
            }
            
            // Create next ray
            ray_i = ray(vertex.rec.p, vertex.wo, ray_i.time());
        }
        
        path.throughput = throughput;
        return path;
    }

    Path generate_light_subpath(const hittable &world, const hittable& lights, int maxDepth)
    {
        Path path;
        color throughput(1, 1, 1);
        
        // Sample a light source
        double lightPdf;
        const hittable& light = lights.sampleHittable(lightPdf);
        
        // Sample a point on the light
        PositionSample positionSample = light.samplePosition(random_double(), random_double());
        
        // Sample emission direction
        DirectionSample directionSample = light.sampleDirection(
            positionSample.p,
            positionSample.n,
            random_double(),
            random_double()
        );
        
        // Create first vertex on light source
        Vertex lightVertex;
        lightVertex.rec.p = positionSample.p;
        lightVertex.rec.normal = positionSample.n;
        lightVertex.rec.mat = light.getMaterial();
        lightVertex.rec.t = 0;
        lightVertex.rec.u = positionSample.u;
        lightVertex.rec.v = positionSample.v;
        // THIS IS A HACK PROBABLY NEED TO THINK OF A SMARTER WAY TO ENFORCE WHICH SIDE IS FRONT
        lightVertex.rec.front_face = true;
        
        lightVertex.wo = directionSample.dir;
        lightVertex.isDelta = false;
        lightVertex.pdfFwd = lightPdf * positionSample.pdf * directionSample.pdf;
        lightVertex.pdfRev = 0.0;
        
        color emission = light.Emission(directionSample.dir, lightVertex.rec);
        double cosTheta = abs(dot(directionSample.dir, positionSample.n));
        throughput = emission * cosTheta / lightVertex.pdfFwd;
        
        lightVertex.throughput = throughput;
        path.vertices.push_back(lightVertex);
        

        ray ray_i(lightVertex.rec.p, directionSample.dir, random_double());
        
        for (int depth = 1; depth < maxDepth; ++depth)
        {
            Vertex vertex;
            
            if (!world.hit(ray_i, interval(0.001, std::numeric_limits<double>::infinity()), vertex.rec))
                break;
            
            vertex.wi = unit_vector(-ray_i.direction());
            
            // Russian roulette
            if (depth > 3)
            {
                double continueProbability = fmin(0.95, color_to_luminance(throughput));
                if (random_double() > continueProbability)
                    break;
                throughput /= continueProbability;
            }
            
            vertex.throughput = throughput;
            
            // Sample BSDF for next direction
            BSDFSample bsdfSample = vertex.rec.mat->Sample(
                vertex.wi,
                vertex.rec,
                random_double(),
                random_double(),
                random_double()
            );
            
            if (bsdfSample.pdf == 0)
                break;
            

            vertex.wo = bsdfSample.dir;
            vertex.isDelta = bsdfSample.isDelta;
            vertex.pdfFwd = bsdfSample.pdf;
            vertex.pdfRev = vertex.rec.mat->Pdf(vertex.wo, vertex.wi, vertex.rec);
            
            path.vertices.push_back(vertex);
            
            if (false)
            {
                throughput = throughput * bsdfSample.val / bsdfSample.pdf;
            }
            else
            {
                cosTheta = abs(dot(vertex.wo, vertex.rec.normal));
                throughput = throughput * bsdfSample.val * cosTheta / bsdfSample.pdf;
            }
            
            ray_i = ray(vertex.rec.p, vertex.wo, ray_i.time());
        }
        
        path.throughput = throughput;
        return path;
    }

    color ConnectPaths(const hittable& world, const hittable& lights, Path cameraPath, Path lightPath)
    {
        // s = light path vertices, t = camera path vertices
        // Try all combinations: 
        // s from 0 to lightPath.vertices.size(), 
        // t from 1 to cameraPath.vertices.size()

        color totalContribution = color(0, 0, 0);
    
        // Optimization: I deleted the MIS Function and am now precomputing all the pdfs
        std::vector<std::vector<double>> pdfs(lightPath.vertices.size() + 1);
        std::vector<std::vector<color>> contributions(lightPath.vertices.size() + 1);
        double sumPdfs = 0.0;
        
        for (int s = 0; s <= lightPath.vertices.size(); ++s)
        {
            pdfs[s].resize(cameraPath.vertices.size() + 1);
            contributions[s].resize(cameraPath.vertices.size() + 1);
            for (int t = 1; t <= cameraPath.vertices.size(); ++t)
            {
                // Calculate contribution with s light vertices and t camera vertices
                contributions[s][t] = EvaluatePath(world, lights, cameraPath, lightPath, s, t);
                if (contributions[s][t].near_zero())
                {
                    pdfs[s][t] = 0.0;
                    continue;
                }
                pdfs[s][t] = PathPdf(world, cameraPath, lightPath, s, t);
                sumPdfs += pdfs[s][t];
            }
        }
        
        if (sumPdfs == 0.0)
        {
            return color(0, 0, 0);
        }

        for (int s = 0; s <= lightPath.vertices.size(); ++s)
        {
           for (int t = 1; t <= cameraPath.vertices.size(); ++t)
           {

               double weight = pdfs[s][t] / sumPdfs;
               // Calculate MIS weight using balance heuristic
               totalContribution += contributions[s][t]*weight;
           }
        }

        return totalContribution;
    }

    color EvaluatePath(const hittable& world, const hittable& lights, Path& cameraPath, Path& lightPath, int s, int t)
    {
        // s = number of light path vertices to use
        // t = number of camera path vertices to use
        
        if (s == 0 && t == 1)
           // Camera directly hits light source
           return EvaluateDirectLightHit(cameraPath);

        if (s == 0 && t > 1)
           // Light path is empty, connect last camera vertex to light
           return ConnectToLight(world, lights, cameraPath, t);
        
        if (s > 0 && t >= 1)
            // Connect camera path vertex t-1 to light path vertex s-1
            return ConnectVertices(world, cameraPath, lightPath, s, t);
        
        return color(0, 0, 0);
    }

    color ConnectVertices(const hittable& world, Path& cameraPath, Path& lightPath, int s, int t)
    {
        Vertex cameraVertex = cameraPath.vertices[t - 1];
        Vertex lightVertex = lightPath.vertices[s - 1];
        
        // Check if connection is possible
        vec3 connectionDir = lightVertex.rec.p - cameraVertex.rec.p;
        double distance = connectionDir.length();
        connectionDir = unit_vector(connectionDir);

        
        // Visibility test
        hit_record rec;
        ray visibility(cameraVertex.rec.p,connectionDir);
        if (world.hit(visibility, interval(0.001, distance-.001), rec))
        {
            return color(0, 0, 0);
        }
        
        // Cannot connect through delta interactions (specular surfaces)
        // if (cameraVertex.isDelta || lightVertex.isDelta)
        // {
        //     return color(0, 0, 0);
        // }
        
        // Evaluate BSDFs at both vertices
        color cameraBSDF = cameraVertex.rec.mat->Evaluate(
            cameraVertex.wi,
            connectionDir,
            cameraVertex.rec
        );
        
        color lightBSDF;
        if (s == 1)
        {
            // Documentation describes this term as the relative emission to 
            // the maximum emission from the light
            // lightBSDF = lightVertex.rec.mat->Emission(
            //     lightVertex.rec.p,
            //     -connectionDir,
            //     lightVertex.rec.normal,
            //     lightVertex.rec.u,
            //     lightVertex.rec.v
            // );
            //For now assume 100% transmission
            lightBSDF = color(1,1,1);
        }
        else
        {
            // Not first light vertex: normal BSDF evaluation
            lightBSDF = lightVertex.rec.mat->Evaluate(
                lightVertex.wi,
                -connectionDir,
                lightVertex.rec
            );
        }

        double cosCameraAngle = abs(dot(connectionDir, cameraVertex.rec.normal));
        double cosLightAngle = abs(dot(-connectionDir, lightVertex.rec.normal));
        
        double geometryTerm = (cosCameraAngle * cosLightAngle) / (distance * distance);
        color contribution = cameraVertex.throughput * cameraBSDF * geometryTerm * 
                lightBSDF * lightVertex.throughput;
        
        return contribution;
    }

    double PathPdf(const hittable& world, Path& cameraPath, Path& lightPath, int s, int t)
    {
        if (t == 0 || 
            t > cameraPath.vertices.size() || 
            s > lightPath.vertices.size())
        {
            return 0.0;
        }
        
        double pdf = 1.0;
        
        // Camera subpath contribution
        if (t >= 1)
        {
            // Camera sampling PDF (stored in first vertex)
            if (cameraPath.vertices[0].pdfRev == 0.0)
                return 0.0;
            pdf *= cameraPath.vertices[0].pdfRev;
            
            // BSDF sampling PDFs
            for (int i = 0; i < t - 1; ++i)
            {
                if (cameraPath.vertices[i].pdfFwd == 0.0)
                    return 0.0;
                pdf *= cameraPath.vertices[i].pdfFwd;
            }
        }
        
        // Light subpath contribution
        if (s >= 1)
        {
            if (lightPath.vertices[0].pdfFwd == 0.0)
                return 0.0;
            pdf *= lightPath.vertices[0].pdfFwd;
            
            // BSDF sampling PDFs
            for (int i = 0; i < s - 1; ++i)
            {
                if (lightPath.vertices[i].pdfFwd == 0.0)
                    return 0.0;
                pdf *= lightPath.vertices[i].pdfFwd;
            }
        }
        
        return pdf;
    }


    color EvaluateDirectLightHit(Path& cameraPath)
    {
        // Camera ray directly hits light source
        if (cameraPath.vertices.size() == 0)
        {
            return color(0, 0, 0);
        }

        // This should be handled by checking if first ray hits emissive surface
        Vertex vertex = cameraPath.vertices[0];
        if (vertex.rec.mat->isEmissive())
        {
            return vertex.rec.mat->Emission(
                vertex.wi,
                vertex.rec);
        }

        return color(0, 0, 0);
    }

    color ConnectToLight(const hittable& world, const hittable& lights, Path& cameraPath, int t)
    {
        Vertex vertex = cameraPath.vertices[t - 1];
        
        // Sample a random light
        double lightPdf;
        const hittable& light = lights.sampleHittable(lightPdf);

        // Sample point on light
        PositionSample positionSample = light.samplePosition(random_double(), random_double());

        // Check visibility and evaluate contribution
        vec3 lightDir = positionSample.p - vertex.rec.p;
        double distance = lightDir.length();
        lightDir = unit_vector(lightDir);


        // Visibility test
        hit_record rec;
        ray visibility(vertex.rec.p,lightDir);
        if (world.hit(visibility, interval(0.001, distance-.001), rec))
        {
            return color(0, 0, 0);
        }
        
        color bsdf = vertex.rec.mat->Evaluate(
            vertex.wi,
            lightDir,
            vertex.rec
        );

        hit_record lrec;
        lrec.p = positionSample.p;
        lrec.normal = positionSample.n;
        lrec.mat = light.getMaterial();
        lrec.t = 0;
        lrec.u = positionSample.u;
        lrec.v = positionSample.v;
        // THIS IS A HACK PROBABLY NEED TO THINK OF A SMARTER WAY TO ENFORCE WHICH SIDE IS FRONT
        lrec.front_face = true;
        

        color emission = light.Emission(
            -lightDir,
            lrec
        );
        
        double cosVertex = abs(dot(lightDir, vertex.rec.normal));
        double cosLight = abs(dot(-lightDir, positionSample.n));
        
        double geometryTerm = (cosVertex * cosLight) / (distance * distance);
        
        return vertex.throughput * bsdf * geometryTerm * emission / 
            (lightPdf * positionSample.pdf);
    }
};

#endif
