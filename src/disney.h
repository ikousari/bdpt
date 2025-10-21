#ifndef DISNEY_H
#define DISNEY_H

#include "material.h"

class disneyClearcoat : public material
{
public:
  disneyClearcoat(double clearcoatGloss) : m_clearcoatGloss{clearcoatGloss} {}

  virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec) const override
  {
    return false;
  }

  virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const override
  {
    return 0;
  }

private:
  double m_clearcoatGloss;
};

class disneyDiffuse : public material
{
public:
  disneyDiffuse(const color &albedo, double roughness, double subsurface)
      : m_baseColor(make_shared<solid_color>(albedo), m_roughness{roughness}, m_subsurface{subsurface}) {}

  disneyDiffuse(shared_ptr<texture> tex, double roughness, double subsurface)
      : m_baseColor(tex), m_roughness{roughness}, m_subsurface{subsurface} {}

  virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec) const override
  {
    return false;
  }

  virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const override
  {
    return 0;
  }

private:
  shared_ptr<texture> m_baseColor;
  double m_roughness;
  double m_subsurface;
};

class disneyGlass : public material
{
public:
  disneyGlass(const color &albedo, double roughness, double anisotropic, double eta)
      : m_baseColor(make_shared<solid_color>(albedo)), m_roughness{roughness}, m_anisotropic{anisotropic}, 
      m_eta{eta} {}

  disneyGlass(shared_ptr<texture> tex, double roughness, double anisotropic, double eta)
      : m_baseColor(tex), m_roughness{roughness}, m_anisotropic{anisotropic}, m_eta{eta} {}

  virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec) const override
  {
    return false;
  }

  virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const override
  {
    return 0;
  }

private:
  shared_ptr<texture> m_baseColor;
  double m_roughness;
  double m_anisotropic;
  double m_eta;
};

class disneyMetal : public material
{
public:
  disneyMetal(const color &albedo, double roughness, double anisotropic)
      : m_baseColor(make_shared<solid_color>(albedo)), m_roughness{roughness}, m_anisotropic{anisotropic} {}

  disneyMetal(shared_ptr<texture> tex, double roughness, double anisotropic)
      : m_baseColor(tex), m_roughness{roughness}, m_anisotropic{anisotropic} {}

  virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec) const override
  {
    return false;
  }

  virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const override
  {
    return 0;
  }

private:
  shared_ptr<texture> m_baseColor;
  double m_roughness;
  double m_anisotropic;
};

class disneySheen : public material
{
public:
  disneySheen(const color &albedo, double sheenTint)
      : m_baseColor(make_shared<solid_color>(albedo)), m_sheenTint{sheenTint} {}

  disneySheen(shared_ptr<texture> tex, double sheenTint)
      : m_baseColor(tex), m_sheenTint{sheenTint} {}

  virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec) const override
  {
    return false;
  }

  virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const override
  {
    return 0;
  }

private:
  shared_ptr<texture> m_baseColor;
  double m_sheenTint;
};

class disney : public material
{
public:
  disney(const color &albedo, double specularTransmission, double metallic,
         double subsurface, double specular, double roughness,
         double specularTint, double anisotropic, double sheen,
         double sheenTint, double clearCoat, double clearcoatGloss, double eta)
  : m_baseColor(make_shared<solid_color>(albedo)), m_specularTransmission{specularTransmission}, 
    m_metallic{metallic}, m_subsurface{subsurface}, m_specular{specular}, m_roughness{roughness},
    m_specularTint{specularTint}, m_anisotropic{anisotropic}, m_sheen{sheen}, m_sheenTint{sheenTint},
    m_clearcoat{clearCoat}, m_clearcoatGloss{clearcoatGloss}, m_eta{eta} {}

  disney(shared_ptr<texture> tex, double specularTransmission, double metallic,
         double subsurface, double specular, double roughness,
         double specularTint, double anisotropic, double sheen,
         double sheenTint, double clearCoat, double clearcoatGloss, double eta)
  : m_baseColor(tex), m_specularTransmission{specularTransmission}, m_metallic{metallic}, 
    m_subsurface{subsurface}, m_specular{specular}, m_roughness{roughness}, 
    m_specularTint{specularTint}, m_anisotropic{anisotropic}, m_sheen{sheen}, 
    m_sheenTint{sheenTint}, m_clearcoat{clearCoat}, m_clearcoatGloss{clearcoatGloss}, m_eta{eta} {}

  virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec) const override
  {
    return false;
  }

  virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const override
  {
    return 0;
  }

private:
  shared_ptr<texture> m_baseColor;
  double m_specularTransmission;
  double m_metallic;
  double m_subsurface;
  double m_specular;
  double m_roughness;
  double m_specularTint;
  double m_anisotropic;
  double m_sheen;
  double m_sheenTint;
  double m_clearcoat;
  double m_clearcoatGloss;
  double m_eta;
};

#endif
