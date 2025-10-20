//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    // 单独设置EPSILON, 避免EPSILON过小导致渲染结果出现条纹
    const float EPSILON = 0.0001f;
    // inter 为光线与物体的交点
    Intersection inter = intersect(ray);
    // inter_material 为物体的材质
    Material* inter_material = inter.m;

    // 如果没有交点
    if (!inter.happened) 
        return Vector3f(0.0f);

    // 如果是光源,返回光源颜色
    if (inter_material->hasEmission()) 
        return inter_material->getEmission();
    
    // 开始计算直接光照
    Vector3f l_dir(0.0f);

    // 这里 wo 为物体到相机(或上一个物体交点)的方向
    // 通常用 wo (outgoing direction)表示交点到相机的方向(出射方向),用 wi (incoming direction)表示交点到光源的方向(入射方向)
    Vector3f wo = -ray.direction;

    // light_inter 为全部光源上的某一随机采样点
    Intersection light_inter;
    // pdf_light 为光源采样点的概率密度函数的值
    float pdf_light = 0.0f;
    sampleLight(light_inter, pdf_light);

    // wi_light 为物体交点到光源的方向
    Vector3f wi_light = light_inter.coords - inter.coords;
    // wi_distance 为物体到光源的距离
    float wi_distance = (light_inter.coords - inter.coords).norm();
    wi_light = normalize(wi_light);
    // light_inter2 为物体往光源方向的交点,用于判断物体和光源间是否有遮挡
    Intersection light_inter2 = intersect(Ray(inter.coords, wi_light));
    // 判断和光源间是否有遮挡 
    if(std::fabs(wi_distance - light_inter2.distance) < EPSILON)
    {
        // f_r 为材质的反射率
        Vector3f f_r = inter_material->eval(wo, wi_light, inter.normal);
        // cos_theta 为交点到光源方向与法线的夹角余弦
        float cos_theta = std::max(0.0f, dotProduct(inter.normal, wi_light));
        // cos_theta_light 为光源到交点方向与光源法线的夹角余弦
        float cos_theta_light = std::max(0.0f, dotProduct(light_inter.normal, -wi_light));
        float distance2 = wi_distance * wi_distance;
        // light_inter.emit 为光源的辐射度
        l_dir = light_inter.emit * f_r * cos_theta * cos_theta_light / distance2 / pdf_light;
    }

    // 开始计算间接光照
    Vector3f l_indir(0.0f);
    // get_random_float() 生成[0,1)之间的随机数
    // 通过俄罗斯轮盘赌决定是否继续追踪间接光照
    if(get_random_float() < RussianRoulette)
    {
        // 根据材质和 wo、法线 N，采样一个新的方向 wi_next，用于间接光照的递归追踪
        Vector3f wi_next = normalize(inter_material->sample(wo, inter.normal));
        Ray next_ray(inter.coords, wi_next);
        Intersection inter_next = intersect(next_ray);
        // 如果有交点且不是光源
        if(inter_next.happened && !inter_next.m->hasEmission())
        {
            // pdf 为材质采样的概率密度函数
            float pdf = inter_material->pdf(wo, wi_next, inter.normal);
            // 如果概率密度函数大于0
            // 避免除以0返回无穷大,导致画面出现亮点
            if(pdf > EPSILON)
            {
                // f_r 为材质的反射率
                Vector3f f_r = inter_material->eval(wo, wi_next, inter.normal);
                float cos_theta = std::max(0.0f, dotProduct(inter.normal, wi_next));
                // 递归计算间接光照
                Ray next_ray(inter.coords, wi_next);
                // castRay(next_ray, depth + 1) 为从交点沿 wi_next 方向追踪的光线的辐射度
                l_indir = castRay(next_ray, depth + 1) * f_r * cos_theta / pdf / RussianRoulette;
            }
        }
 
    }

    // 返回直接光照和间接光照之和
    return l_dir + l_indir;
}   