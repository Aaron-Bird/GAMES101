#include <algorithm>
#include <cassert>
#include "BVH.hpp"
#include <chrono>

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    auto t_start = std::chrono::high_resolution_clock::now();
    root = recursiveBuild(primitives);
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    printf("BVH build time: %.3f ms\n", elapsed);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        std::vector<Object*> leftshapes;
        std::vector<Object*> rightshapes;

        if(splitMethod == SplitMethod::SAH && objects.size() > 4)
        {
            for (int i = 0; i < objects.size(); i++)
                centroidBounds = Union(centroidBounds, objects[i]->getBounds().Centroid());

            int bestAxis = 0;
            int bestSplit = 0;
            float minCost = std::numeric_limits<float>::infinity();
            const int numBuckets = 12;

            std::vector<std::array<int, 3>> objectBucketIndices(objects.size());
            for (int i = 0; i < objects.size(); i++) 
            {
                std::array<int, 3> bucketIndex;
                Vector3f offset = centroidBounds.Offset(objects[i]->getBounds().Centroid());
                bucketIndex[0] = std::clamp(int(numBuckets * offset.x), 0, numBuckets - 1);
                bucketIndex[1] = std::clamp(int(numBuckets * offset.y), 0, numBuckets - 1);
                bucketIndex[2] = std::clamp(int(numBuckets * offset.z), 0, numBuckets - 1);
                objectBucketIndices[i] = bucketIndex;
            }

            for( int axis = 0; axis < 3; axis++) 
            {
                int bucketCounts[numBuckets] {0};
                Bounds3 bucketBounds[numBuckets];
               
                for (int i = 0; i < objects.size(); i++) 
                {
                    // int bucketIndex = numBuckets * centroidBounds.Offset(objects[i]->getBounds().Centroid())[axis];
                    // bucketIndex = std::clamp(bucketIndex, 0, numBuckets - 1);
                    int bucketIndex = objectBucketIndices[i][axis];
                    bucketCounts[bucketIndex]++;
                    bucketBounds[bucketIndex] = Union(bucketBounds[bucketIndex], objects[i]->getBounds());
                }

                Bounds3 prefixBounds[numBuckets], suffixBounds[numBuckets];
                int prefixCounts[numBuckets] {0}, suffixCounts[numBuckets] {0};
                prefixBounds[0] = bucketBounds[0];
                prefixCounts[0] = bucketCounts[0];
                suffixBounds[numBuckets - 1] = bucketBounds[numBuckets - 1];
                suffixCounts[numBuckets - 1] = bucketCounts[numBuckets - 1];

                for (int i = 1; i < numBuckets; ++i)
                {
                    prefixBounds[i] = Union(prefixBounds[i-1], bucketBounds[i]);
                    prefixCounts[i] = prefixCounts[i-1] + bucketCounts[i];
                }
                
                for (int i = numBuckets - 2; i >= 0; i--)
                {
                    suffixBounds[i] = Union(suffixBounds[i+1], bucketBounds[i]);
                    suffixCounts[i] = suffixCounts[i+1] + bucketCounts[i];
                }

                for(int i = 0; i < numBuckets - 1; i++) 
                {
                    Bounds3 leftBound = prefixBounds[i];
                    Bounds3 rightBound = suffixBounds[i+1];
                    int leftCount = prefixCounts[i];
                    int rightCount = suffixCounts[i+1];

                    float cost = 1 + (leftCount * leftBound.SurfaceArea() + rightCount * rightBound.SurfaceArea()) / bounds.SurfaceArea();

                    if(cost < minCost) 
                    {
                        minCost = cost;
                        bestAxis = axis;
                        bestSplit = i;
                    }
                }
            }

            for (int i = 0; i < objects.size(); i++) 
            {
                int bucketIndex = objectBucketIndices[i][bestAxis];
                if (bucketIndex <= bestSplit)
                    leftshapes.push_back(objects[i]);
                else
                    rightshapes.push_back(objects[i]);
            }
        } else {
            for (int i = 0; i < objects.size(); ++i)
                centroidBounds = Union(centroidBounds, objects[i]->getBounds().Centroid());
            int dim = centroidBounds.maxExtent();
            switch (dim) {
            case 0:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().x <
                           f2->getBounds().Centroid().x;
                });
                break;
            case 1:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().y <
                           f2->getBounds().Centroid().y;
                });
                break;
            case 2:
                std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                    return f1->getBounds().Centroid().z <
                           f2->getBounds().Centroid().z;
                });
                break;
            }

            auto beginning = objects.begin();
            auto middling = objects.begin() + (objects.size() / 2);
            auto ending = objects.end();

            leftshapes = std::vector<Object*>(beginning, middling);
            rightshapes = std::vector<Object*>(middling, ending);
        }   
     
    

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    std::array<int, 3> dirIsNeg { 
        int(ray.direction.x > 0), 
        int(ray.direction.y > 0), 
        int(ray.direction.z > 0) 
    };
    if(!node || !node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg ))
        return Intersection();

    if(node->left == nullptr && node->right == nullptr && node->object)
        return node->object->getIntersection(ray);

    Intersection leftInter = getIntersection(node->left, ray);
    Intersection rightInter = getIntersection(node->right, ray);

    if(leftInter.happened && rightInter.happened) 
        return leftInter.distance < rightInter.distance ? leftInter : rightInter;
    else if(leftInter.happened)
        return leftInter;
    else if(rightInter.happened)
        return rightInter;
    else
        return Intersection();
}