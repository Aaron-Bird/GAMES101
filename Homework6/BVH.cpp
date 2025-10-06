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
        // centroidBounds 为所有 object 的中心点的包围盒
        Bounds3 centroidBounds;
        // leftshapes 和 rightshapes 分别为划分后左右子节点的 object 列表
        std::vector<Object*> leftshapes;
        std::vector<Object*> rightshapes;
        // 仅当 object 数量 > 4 时使用 SAH 划分算法
        // object 数量 <= 4 时使用中位数划分算法
        if(splitMethod == SplitMethod::SAH && objects.size() > 4)
        {

            for (int i = 0; i < objects.size(); i++)
                centroidBounds = Union(centroidBounds, objects[i]->getBounds().Centroid());
            // bestAxis, bestSplit 分别为最佳划分轴和划分位置
            int bestAxis = 0;
            int bestSplit = 0;
            // minCost 为当前最小划分代价
            float minCost = std::numeric_limits<float>::infinity();
            // 因为对 objects 排序耗时较多,使用桶排序来减少排序时间
            // 根据 object 的中心点位置将其划分到不同的桶中
            // numBuckets 为划分桶的数量
            const int numBuckets = 12;

            // objectBucketIndices 为每个 object 所在的桶的索引
            // 放在 for 循环外面是为了供后续的 leftshapes.push_back 和 rightshapes.push_back 使用
            std::vector<std::array<int, 3>> objectBucketIndices(objects.size());
            for (int i = 0; i < objects.size(); i++) 
            {
                std::array<int, 3> bucketIndex;
                // 根据 object 的中心点到 centroidBounds 的偏移量来计算其所在的桶的索引
                Vector3f offset = centroidBounds.Offset(objects[i]->getBounds().Centroid());
                // 使用 clamp 是为了确保 bucketIndex 在 [0, numBuckets - 1] 范围内
                bucketIndex[0] = std::clamp(int(numBuckets * offset.x), 0, numBuckets - 1);
                bucketIndex[1] = std::clamp(int(numBuckets * offset.y), 0, numBuckets - 1);
                bucketIndex[2] = std::clamp(int(numBuckets * offset.z), 0, numBuckets - 1);
                objectBucketIndices[i] = bucketIndex;
            }

            // 遍历三个轴,分别计算在该轴上划分的代价
            for( int axis = 0; axis < 3; axis++) 
            {
                // bucketCounts 为每个桶中包含的 object 数量
                int bucketCounts[numBuckets] {0};
                // bucketBounds 为每个桶对应的包围盒
                Bounds3 bucketBounds[numBuckets];
               
                for (int i = 0; i < objects.size(); i++) 
                {
                    // int bucketIndex = numBuckets * centroidBounds.Offset(objects[i]->getBounds().Centroid())[axis];
                    // bucketIndex = std::clamp(bucketIndex, 0, numBuckets - 1);
                    // 取出该 object 在该轴上的桶索引
                    int bucketIndex = objectBucketIndices[i][axis];

                    // 更新对应桶的 object 数量和包围盒
                    bucketCounts[bucketIndex]++;
                    bucketBounds[bucketIndex] = Union(bucketBounds[bucketIndex], objects[i]->getBounds());
                }

                // prefixBounds 和 suffixBounds 为当前桶索引 i 的左侧和右侧所有桶的包围盒
                Bounds3 prefixBounds[numBuckets], suffixBounds[numBuckets];
                // prefixCounts 和 suffixCounts 为当前桶索引 i 的左侧和右侧所有桶的 object 数量
                int prefixCounts[numBuckets] {0}, suffixCounts[numBuckets] {0};
                prefixBounds[0] = bucketBounds[0];
                prefixCounts[0] = bucketCounts[0];
                suffixBounds[numBuckets - 1] = bucketBounds[numBuckets - 1];
                suffixCounts[numBuckets - 1] = bucketCounts[numBuckets - 1];
                
                for (int i = 1; i < numBuckets; ++i)
                {
                    // prefixBounds 包含自身和左侧所有桶的包围盒
                    // prefixCounts 包含自身和左侧所有桶的 object 数量
                    prefixBounds[i] = Union(prefixBounds[i-1], bucketBounds[i]);
                    prefixCounts[i] = prefixCounts[i-1] + bucketCounts[i];
                }
                
                for (int i = numBuckets - 2; i >= 0; i--)
                {
                    // suffixBounds 包含自身和右侧所有桶的包围盒
                    // suffixCounts 包含自身和右侧所有桶的 object 数量
                    suffixBounds[i] = Union(suffixBounds[i+1], bucketBounds[i]);
                    suffixCounts[i] = suffixCounts[i+1] + bucketCounts[i];
                }

                // 从左向右遍历桶,将桶作为划分点,计算划分代价
                for(int i = 0; i < numBuckets - 1; i++) 
                {
                    Bounds3 leftBound = prefixBounds[i];
                    Bounds3 rightBound = suffixBounds[i+1];
                    int leftCount = prefixCounts[i];
                    int rightCount = suffixCounts[i+1];
                    // 计算当前划分的代价
                    // 此处假设 t_traversal = 1, t_intersect = 1
                    float cost = 1 + (leftCount * leftBound.SurfaceArea() + rightCount * rightBound.SurfaceArea()) / bounds.SurfaceArea();

                    if(cost < minCost) 
                    {
                        minCost = cost;
                        bestAxis = axis;
                        bestSplit = i;
                    }
                }
            }

            // 根据最佳划分轴和最佳划分位置,将 objects 放入对应的 leftshapes 和 rightshapes
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